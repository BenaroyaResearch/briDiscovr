## Copyright (C) 2020  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the briDiscovr package

#' Prepare a dataset for DISCOV-R analysis
#'
#' Loads in the marker and .fcs information and prepares a \code{discovrExperiment}
#' object for analysis. Default parameters values are based on the original DISCOV-R
#' analysis published in Wiedeman et al 2020.
#'
#' @param markerInfoFile A character string indicating the path to a .csv file. This file is expected to have columns
#' named useToCluster", as well as the names specified in the \code{markerCommonField} and \code{markerFcsField} variables.
#' Details of the \code{markerCommonField} and \code{markerFcsField} arguments are provided below; the "useToCluster"
#' column should have only TRUE or FALSE values. Markers with a TRUE value in this column will be used for clustering,
#' whereas the others will not.
#' @param fcsInfoFile A character string indicating the path to a file containing columns named "subject",
#' "cellSubset", and "filename". The "filename" field must contain paths to the .fcs files that will be used in analysis.
#' @param markerCommonField A character string indicating the name of a column containing common marker names for human use,
#' like "CD45" (default: "fixed")
#' @param markerFcsField A character string indicating the name of a column containing the marker names in the .fcs files,
#' like "89Y_CD45" (default: "desc")
#' @param arcsinhA (default: 0) A numeric indicating the value for 'a' in the arcsinh data transformation
#' equation. Should usually be 0.
#' @param arcsinhB (default: 0.2) A numeric indicating the value for 'b' in the arcsinh data transformation
#' equation. Should be 1/5 for cytof data, and 1/150 for flow data.
#' @param arcsinhC (default: 0) A numeric indicating the value for 'c' in the arcsinh data transformation
#' equation. Should usually be 0.
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{discovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import flowCore
#' @import dplyr
#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
#' @export
setupDiscovrExperiment <- function(
  markerInfoFile,
  fcsInfoFile,
  markerCommonField = "fixed",
  markerFcsField = "desc",
  arcsinhA = 0,
  arcsinhB = 0.2,
  arcsinhC = 0,
  verbose = TRUE
){
  # get marker info, and check for appropriate columns
  markerInfo <- read.csv(markerInfoFile, stringsAsFactors = FALSE)
  if(!all(c(markerCommonField, markerFcsField, "useToCluster") %in% names(markerInfo))){
    stop(
      "The file set as 'markerInfoFile' must contain columns with names 'useToCluster' and the names
      specified by commonMarkerField and fcsMarkerField. Please check your markerInfo file and try again."
    )
  }
  if(!all(is.logical(markerInfo$useToCluster))){
    stop("The 'useToCluster' column in the 'markerInfo' file must contain only TRUE or FALSE values.")
  }

  clusteringMarkers <- markerInfo[markerInfo$useToCluster, "commonMarkerName", drop = TRUE]
  if(verbose){message(paste("Found clustering markers:", paste(clusteringMarkers, collapse=", ")))}

  # get list of FCS files with subject and cell subset information, check for appropriate columns
  fcsInfo = read.csv(fcsInfoFile, stringsAsFactors = FALSE)
  if(!all(c("subject", "cellSubset", "filename") %in% names(fcsInfo))){
    stop("The file set as 'fcsInfoFile' must contain columns with names 'subject', 'cellSubset', and 'filename'.")
  }

  # Check that fcs files can be found
  filesThatDontExist <- c()
  for(currFile in fcsInfo$filename){
    if(!file.exists(currFile)){
      filesThatDontExist <- c(filesThatDontExist, currFile)
    }
  }
  if(length(filesThatDontExist > 0)){
    stop(paste(
      "The following files could not be found:\n",
      paste0(filesThatDontExist, collapse = "\n")
    ))
  }

  # Check fcs files for byte offset issue that will prevent analysis
  filesWithOffsetIssues <- c()
  for(currFile in fcsInfo$filename){
    if(checkForFcsByteOffsetIssue(currFile)){
      filesWithOffsetIssues <- c(filesWithOffsetIssues, currFile)
    }
  }
  if(length(filesWithOffsetIssues > 0)){
    stop(paste(
      "The following files have byte offset issues and cannot be used. Please try re-exporting these files:\n",
      paste0(filesWithOffsetIssues, collapse = "\n")
    ))
  }

  # Check fcs files for 0 events, which will prevent analysis
  filesWithNoEvents <- c()
  for(currFile in fcsInfo$filename){
    if(!getFcsNEvents(currFile) > 0){
      filesWithNoEvents <- c(filesWithNoEvents, currFile)
    }
  }
  if(length(filesWithNoEvents > 0)){
    stop(paste(
      "The following files have no events and cannot be used. Please remove these from your list of .fcs files:\n",
      paste0(filesWithNoEvents, collapse = "\n")
    ))
  }

  if(verbose){message("All .fcs files passed basic QC checks. Proceeding...")}

  ###########################################################
  # Section 2.c.i-iii from original SOP load and format input
  ###########################################################

  # Read in FCS files and associate with cell subsets and subjects
  for (currCellSubset in unique(fcsInfo$cellSubset)){
    if(verbose){message(paste0("Assigning per-subject data for subset: ", currCellSubset, "..."))}
    assign(currCellSubset, buildFcsList(fcsInfo %>% dplyr::filter(cellSubset == currCellSubset), truncate_max_range = FALSE))
  }
  assign("All_Subjects", buildFcsList(fcsInfo[!duplicated(fcsInfo$subject),], truncate_max_range = FALSE))

  # Process data
  processData <- function(fcs){
    # Tidy marker names
    pData(parameters(fcs))$desc <-
      markerInfo$commonMarkerName[match(pData(parameters(fcs))$desc, markerInfo$fcsMarkerName)]

    # This changes parameters(fcs)$name, featureNames(fcs), and colnames(fcs) - aka events colnames - all in one fell swoop.
    # note colnames has to be the one from flowCore
    flowCore::colnames(fcs) = make.names(pData(parameters(fcs))$desc)

    # Remove markers that aren't informative/shared between panels (i.e. duplicated NAs)
    fcs = fcs[,!(duplicated(flowCore::colnames(fcs)) | duplicated(flowCore::colnames(fcs), fromLast = TRUE))]
    fcs = fcs[, order(flowCore::colnames(fcs))]
  }

  if(verbose){message("Cleaning and selecting markers from .fcs files according to markerInfo file...")}

  for (currCellSubset in unique(fcsInfo$cellSubset)){
    tmpFile <- get(currCellSubset)
    tmpFile <- lapply(tmpFile, processData)
    assign(currCellSubset, tmpFile)
  }
  All_Subjects <- lapply(All_Subjects, processData)

  if(verbose){message("Merging data (parent and child populations)...")}
  allMerged = All_Subjects
  for(currSubject in names(All_Subjects)){
    exprs(allMerged[[currSubject]]) = exprs(allMerged[[currSubject]])[0,]
    for(cellSubset in unique(fcsInfo$cellSubset)){
      exprs(allMerged[[currSubject]]) = rbind(
        exprs(allMerged[[currSubject]]),
        if(!is.null(get(cellSubset)[[currSubject]])) exprs(get(cellSubset)[[currSubject]])
      )
    }
  }

  # Transform the data
  asinhTfmData <- function(fcs){
    # Arcsinh transform remaining columns
    tl <- transformList(
      flowCore::colnames(fcs),
      arcsinhTransform(a=arcsinhA, b=arcsinhB, c=arcsinhC),
      transformationId="asinh"
    )
    fcs = transform(fcs, tl)
  }

  if(verbose){
    message(paste0("Transforming data using arcsinh(a=", arcsinhA, ", b=", arcsinhB, ", c=", arcsinhC, ")..."))
  }
  allDataTransformed <- lapply(allMerged, asinhTfmData)
  transformedFlowSet <- as(allDataTransformed, "flowSet")

  # Extract expression data and label Tmr+ events
  mergedExpr <- setNames(
    data.frame(matrix(ncol = ncol(exprs(transformedFlowSet[[1]]))+2, nrow = 0)),
    c(flowCore::colnames(transformedFlowSet), "samp", "cellSubset")
  )

  for(currSubject in names(All_Subjects)){
    tmpExpr = as.data.frame(flowCore::exprs(transformedFlowSet[[currSubject]]))
    tmpExpr$samp = as.character(currSubject)
    temp_cellSubsets <- list()
    for(cellSubset in unique(fcsInfo$cellSubset)){
      temp_cellSubsets <- c(
        temp_cellSubsets,
        if(!is.null(get(cellSubset)[[currSubject]])) rep(paste0(cellSubset), nrow(flowCore::exprs(get(cellSubset)[[currSubject]])))
      )
    }
    tmpExpr$cellSubset = temp_cellSubsets
    mergedExpr = rbind(mergedExpr, tmpExpr)
  }

  # building a discovr experiment S3 object
  exptInProgress <- structure(list(), class = "discovrExperiment")
  exptInProgress$markerInfoFile <- markerInfoFile
  exptInProgress$markerInfo     <- markerInfo
  exptInProgress$fcsInfoFile    <- fcsInfoFile
  exptInProgress$fcsInfo        <- fcsInfo
  exptInProgress$mergedExpr     <- mergedExpr
  return(exptInProgress)
}

#' Load in a dataset and perform phenograph clustering
#'
#' @param markerInfoFile A character string indicating the path to a .csv file. This file is expected to have columns
#' named useToCluster", as well as the names specified in the \code{markerCommonField} and \code{markerFcsField} variables.
#' Details of the \code{markerCommonField} and \code{markerFcsField} arguments are provided below; the "useToCluster"
#' column should have only TRUE or FALSE values. Markers with a TRUE value in this column will be used for clustering,
#' whereas the others will not.
#' @param fcsInfoFile A character string indicating the path to a file containing columns named "subject",
#' "cellSubset", and "filename". The "filename" field must contain paths to the .fcs files that will be used in analysis.
#' @param markerCommonField A character string indicating the name of a column containing common marker names for human use,
#' like "CD45" (default: "fixed")
#' @param markerFcsField A character string indicating the name of a column containing the marker names in the .fcs files,
#' like "89Y_CD45" (default: "desc")
#' @param arcsinhA (default: 0) A numeric indicating the value for 'a' in the arcsinh data transformation
#' equation. Should usually be 0.
#' @param arcsinhB (default: 0.2) A numeric indicating the value for 'b' in the arcsinh data transformation
#' equation. Should be 1/5 for cytof data, and 1/150 for flow data.
#' @param arcsinhC (default: 0) A numeric indicating the value for 'c' in the arcsinh data transformation
#' equation. Should usually be 0.
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import flowCore
#' @import dplyr
#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
#' @export
setupDiscovrExperiment <- function(
  markerInfoFile,
  fcsInfoFile,
  markerCommonField = "fixed",
  markerFcsField = "desc",
  arcsinhA = 0,
  arcsinhB = 0.2,
  arcsinhC = 0,
  verbose = TRUE
){

  ###################################################################
  # 2.d.i: Phenograph clustering - takes a while
  ###################################################################
  # Clustering markers object is defined in the first set-up chunk.  Tweak which markers are included up there.
#
#   # Set up Phenograph function to use kd tree
#   find_neighbors <- function(data, k){
#     nearest <- RANN::nn2(data, data, k, treetype = "kd", searchtype = "standard")
#     return(nearest[[1]])
#   }
#
#   Rpheno <- function(data, k=30){
#     if(is.data.frame(data))
#       data <- as.matrix(data)
#
#     if(!is.matrix(data))
#       stop("Wrong input data, should be a data frame or matrix!")
#
#     if(k<1){
#       stop("k must be a positive integer!")
#     }else if (k > nrow(data)-2){
#       stop("k must be smaller than the total number of points!")
#     }
#
#     message("Run Rphenograph starts:","\n",
#       "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
#       "  -k is set to ", k)
#
#     cat("  Finding nearest neighbors...")
#     t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
#     cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
#     t2 <- system.time(links <- Rphenograph:::jaccard_coeff(neighborMatrix))
#
#     cat("DONE ~",t2[3],"s\n", " Build undirected graph from the weighted links...")
#     links <- links[links[,1]>0, ]
#     relations <- as.data.frame(links)
#     colnames(relations)<- c("from","to","weight")
#     t3 <- system.time(g <- igraph::graph.data.frame(relations, directed=FALSE))
#
#     cat("DONE ~",t3[3],"s\n", " Run louvain clustering on the graph ...")
#     t4 <- system.time(community <- igraph::cluster_louvain(g))
#     cat("DONE ~",t4[3],"s\n")
#
#     message("Run Rphenograph DONE, took a total of ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
#     cat("  Return a community class\n  -Modularity value:", igraph::modularity(community),"\n")
#     cat("  -Number of clusters:", length(unique(igraph::membership(community))))
#
#     return(community)
#   }
#
#   # Run phenograph (using kd treetype) on each subject.
#   PhenographClust = function(fcs, clusteringMarkers) {
#     exprs_mat = as.matrix(as.data.frame(flowCore::exprs(fcs))[,clusteringMarkers])
#     RPvect = as.numeric(igraph::membership(Rpheno(data = exprs_mat)))
#     return(RPvect)
#   }
#   mergedExpr$RPclust = unlist(lapply(allDataTransformed, PhenographClust, clusteringMarkers))
#
#   # Get summary of the number of clusters generated for each subject
#   n_pheno_clusts <- mergedExpr %>%
#     group_by(samp) %>%
#     summarize(k_clusters = max(RPclust))
#   View(n_pheno_clusts)
#
#   ###################################################################
#   # 2.e.i: Summarize and save outputs
#   ###################################################################
#   # Calculate mean expression value of each marker for each phenograph cluster in each subject
#   RP_mean <- mergedExpr %>%
#     dplyr::select(-cellSubset, -contains("tsne"), -contains("umap")) %>%
#     group_by(samp, RPclust) %>%
#     summarise_all(mean) %>%
#     mutate(RPclust = as.character(RPclust))
#
#   # Calculate total CD8 mean expression for each subject  ##NEED TO CHANGE -cell_subset to somethinge else???
#   CD8_mean = mergedExpr %>%
#     dplyr::select(-cellSubset, -RPclust, -contains("tsne"), -contains("umap")) %>%
#     group_by(samp) %>%
#     summarise_all(mean) %>%
#     mutate(RPclust = "Total_CD8")
#
#   RP_mean = bind_rows(RP_mean, CD8_mean)
#
#   # Count cells of each specificity in each phenograph cluster (from each sample)
#
#   RPtmr_counting = mergedExpr %>%
#     dplyr::select(samp, cellSubset, RPclust) %>%
#     group_by(samp, RPclust) %>%
#     summarise(Total=n())
#   RPtmr_counting=as.data.frame(RPtmr_counting)
#
#   uniqueCellSubsets=unique(mergedExpr$cellSubset)
#
#   for(currCellSubset in unique(mergedExpr$cellSubset)){
#     additional_matrix = mergedExpr %>%
#       dplyr::select(samp, cellSubset, RPclust) %>%
#       group_by(samp, RPclust) %>%
#       summarise(number = sum(cellSubset == currCellSubset))
#     additional_matrix=as.data.frame(additional_matrix)
#     print(cellSubset)
#     RPtmr_counting<-cbind(RPtmr_counting, additional_matrix$number)
#   }
#
#   for(i in 1:(ncol(RPtmr_counting)-3)){
#     colnames(RPtmr_counting)[i+3]= (unique(mergedExpr$cellSubset))[i]
#   }
#
#
#   aggregate_counts = RPtmr_counting %>%
#     dplyr::select(-RPclust) %>%
#     group_by(samp) %>%
#     summarise_all(sum) %>%
#     rename_at(vars(-samp),function(name) paste0(name,"_tot"))
#
#   RPtmr_counting = RPtmr_counting %>%
#     left_join(aggregate_counts)
#
#
#   for(i in 1:(length(unique(mergedExpr$cellSubset))+1)){
#     RPtmr_counting = cbind(RPtmr_counting,  RPtmr_counting[,i+2]/RPtmr_counting[,i+2+length(unique(mergedExpr$cellSubset))+1]*100)
#   }
#
#   for(i in 1:(length(unique(mergedExpr$cellSubset)))){
#     colnames(RPtmr_counting)[i+3+(length(unique(mergedExpr$cellSubset))+1)*2]= (paste0("pct_", (unique(mergedExpr$cellSubset))[i], "_in_clust"))
#   }
#
#   colnames(RPtmr_counting)[3+(length(unique(mergedExpr$cellSubset))+1)*2]= "pct_Total_in_clust"
#
#   # Generate file structure in which to save R data and exports
#   outputPathBase <- file.path(workingDir, "Output")
#   dir.create(outputPathBase)
#   dir.create(file.path(outputPathBase, format(Sys.Date(), "%y%m%d")))
#   dir.create(file.path(outputPathBase, format(Sys.Date(), "%y%m%d"), "R_objects"))
#
#   fnamePrefixR <- file.path(outputPathBase, format(Sys.Date(), "%y%m%d"), "R_objects", format(Sys.Date(), "%y%m%d"))
#
#   # Save the clustering output
#   pheno_filename = paste0(fnamePrefixR, "_all_phenograph_data.RData")
#   save(mergedExpr, clusteringMarkers, RP_mean, RPtmr_counting,
#     file = pheno_filename)
#
# }
