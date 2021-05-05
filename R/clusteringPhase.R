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
#' @param parentPopulation A character sting indicating the name of the parent population subset. Must match
#' one of the values in the fcsInfoFile 'cellSubset' column.
#' @param markerCommonField (default: "fixed") A character string indicating the
#' name of a column containing common marker names for human use, like "CD45"
#' @param markerFcsField (default: "desc") A character string indicating the
#' name ofa column containing the marker names in the .fcs files,like "89Y_CD45"
#' @param arcsinhA (default: 0) A numeric indicating the value for 'a' in the
#' arcsinh data transformation equation. Should usually be 0.
#' @param arcsinhB (default: 0.2) A numeric indicating the value for 'b' in the arcsinh
#' data transformation equation. Should be 1/5 for cytof data, and 1/150 for flow data.
#' @param arcsinhC (default: 0) A numeric indicating the value for 'c' in the
#' arcsinh data transformation equation. Should usually be 0.
#' @param verbose (default: TRUE) A boolean specifying whether to display processing messages
#' @param checkMemory (default: TRUE) A boolean indicating whether to check how much system memory
#' is available before loading the dataset. If TRUE, this function will display a message and prevent
#' data loading when the files take up more than 80 percent of the available system memory.
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{discovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import flowCore
#' @import dplyr
#' @importFrom methods as
#' @importFrom stats setNames
#' @importFrom utils read.csv
#' @importFrom memuse Sys.meminfo mu.size
#' @export
setupDiscovrExperiment <- function(
  markerInfoFile,
  fcsInfoFile,
  parentPopulation,
  markerCommonField = "fixed",
  markerFcsField = "desc",
  arcsinhA = 0,
  arcsinhB = 0.2,
  arcsinhC = 0,
  verbose = TRUE,
  checkMemory = TRUE
){
  # get marker info, and check for appropriate columns
  markerInfo <- read.csv(markerInfoFile, stringsAsFactors = FALSE)
  if(!all(c(markerCommonField, markerFcsField, "useToCluster") %in% names(markerInfo))){
    stop(
      "The file set as 'markerInfoFile' must contain columns with names 'useToCluster' and the names ",
      "specified by markerCommonField and markerFcsField Please check your markerInfo file and try again."
    )
  }
  if(!all(is.logical(markerInfo$useToCluster))){
    stop("The 'useToCluster' column in the 'markerInfo' file must contain only TRUE or FALSE values.")
  }
  # set colnames of marker info data for consistency
  markerInfo <- dplyr::rename(markerInfo, commonMarkerName = !!markerCommonField, fcsMarkerName = !!markerFcsField)

  # check for appropriate marker names
  if(nrow(markerInfo) > length(unique(markerInfo$commonMarkerName))){
    stop("Marker names in the ", markerCommonField, " column of the markerInfo file do not appear to be unique. Please correct this and try again.")
  }
  if(any(markerInfo$commonMarkerName == "")){
    stop("Markers without assigned names were detected in the markerInfo file. Please correct the file and try again.")
  }

  # extract clustering markers for convenience
  clusteringMarkers <- markerInfo[markerInfo$useToCluster, "commonMarkerName", drop = TRUE]
  if(verbose){message(paste("Found clustering markers:", paste(clusteringMarkers, collapse=", ")))}

  # get list of FCS files with subject and cell subset information, check for appropriate columns
  fcsInfo = read.csv(fcsInfoFile, stringsAsFactors = FALSE)
  if(!all(c("subject", "cellSubset", "filename") %in% names(fcsInfo))){
    stop("The file set as 'fcsInfoFile' must contain columns with names 'subject', 'cellSubset', and 'filename'.")
  }

  # check to ensure marker names and cell subsets are non-overlapping. Creates issues downstream if this happens.
  overlapNames <- fcsInfo$cellSubset[fcsInfo$cellSubset %in% markerInfo$commonMarkerName]
  if (length(overlapNames) > 0){
    stop(
      "Cell subset names and marker names must not overlap. The following names were detected in both sets of names: ",
      paste0(overlapNames, collapse = ", ")
    )
  }

  # Check that fcs files can be found, and get size of all files
  filesThatDontExist <- c()
  sumFileSize <- 0
  for(currFile in fcsInfo$filename){
    if(!file.exists(currFile)){
      filesThatDontExist <- c(filesThatDontExist, currFile)
    } else {
      sumFileSize <- sumFileSize + file.info(currFile)$size
    }
  }
  if(length(filesThatDontExist > 0)){
    stop(
      "The following files could not be found:\n",
      paste0(filesThatDontExist, collapse = "\n")
    )
  }

  # Check that the files don't exceed the available system memory
  currMem <- memuse::Sys.meminfo()
  currSysMem <- mu.size(currMem$totalram, as.is = FALSE)
  if(sumFileSize > currSysMem){
    stop(
      "Total file size of dataset (", sumFileSize,
      " bytes) exceeds the total system memory (", currSysMem, " bytes)."
    )
  }
  currFreeMem <- mu.size(currMem$freeram, as.is = FALSE)
  if(checkMemory & (sumFileSize > 0.8*currFreeMem)){
    stop(
      "Total file size of dataset (", sumFileSize,
      " bytes) exceeds 80% of the available memory (", currFreeMem, " bytes). ",
      "To process the dataset anyway, you can set the 'checkMemory' option to FALSE, ",
      "but be aware this may cause issues on your system if you run out of memory."
    )
  }
  if(verbose){
    message("Total system memory: ", currSysMem, " bytes")
    message("Available system memory: ", currFreeMem, " bytes")
    message("Size of dataset on disk: ", sumFileSize, " bytes - OK")
  }

  # check that the parent population is one of the subsets in the .fcs info
  if(!parentPopulation %in% unique(fcsInfo$cellSubset)){
    stop(
      "The requested parent population '",
      parentPopulation,
      "' was not found in the 'cellSubset' column of the .fcs info file."
    )
  }

  # check that each subject has only one .fcs file for each cellSubset
  duplicatedSamples = duplicated(paste0(fcsInfo$subject, fcsInfo$cellSubset))
  if(any(duplicatedSamples)){
    stop(
      "The following subject/cell subset pairs were duplicated in the .fcs list:\n",
      apply(
        fcsInfo[duplicatedSamples, c("subject", "cellSubset")],
        1,
        function(r){paste(c(r[1], " - ", r[2], "\n"))}
      )
    )
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
    assign(currCellSubset, buildFcsList(fcsInfo %>% dplyr::filter(.data$cellSubset == currCellSubset), truncate_max_range = FALSE))
  }
  assign("allSubjects", buildFcsList(fcsInfo[fcsInfo$cellSubset == parentPopulation,], truncate_max_range = FALSE))

  # Process data
  processData <- function(fcs){
    # confirm that all markers to be used in clustering are mappable to fcs parameter names
    clusteringMarkerDesc <- markerInfo[markerInfo$useToCluster, "fcsMarkerName", drop = TRUE]
    missingMarkers <- clusteringMarkerDesc[!clusteringMarkerDesc %in% pData(parameters(fcs))$desc]
    if(length(missingMarkers) > 0){
      stop("Could not find the following clustering markers in the .fcs data\n", paste0(missingMarkers, collapse = ", "))
    }

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
  allSubjects <- lapply(allSubjects, processData)

  if(verbose){message("Merging data (parent and child populations)...")}
  allMerged = allSubjects
  for(currSubject in names(allSubjects)){
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
      arcsinhTransform(a = arcsinhA, b = arcsinhB, c = arcsinhC),
      transformationId="asinh"
    )
    fcs = transform(fcs, tl)
  }

  if(verbose){
    message("Transforming data using arcsinh(a=", arcsinhA, ", b=", arcsinhB, ", c=", arcsinhC, ")...")
  }
  allDataTransformed <- lapply(allMerged, asinhTfmData)
  transformedFlowSet <- as(allDataTransformed, "flowSet")

  # Extract expression data and label Tmr+ events
  mergedExpr <- setNames(
    data.frame(matrix(ncol = ncol(exprs(transformedFlowSet[[1]]))+2, nrow = 0)),
    c(flowCore::colnames(transformedFlowSet), "samp", "cellSubset")
  )

  for(currSubject in names(allSubjects)){
    tmpExpr = as.data.frame(flowCore::exprs(transformedFlowSet[[currSubject]]))
    tmpExpr$samp = as.character(currSubject)
    temp_cellSubsets <- list()
    for(cellSubset in unique(fcsInfo$cellSubset)){
      temp_cellSubsets <- c(
        temp_cellSubsets,
        if(!is.null(get(cellSubset)[[currSubject]])){
          rep(paste0(cellSubset), nrow(flowCore::exprs(get(cellSubset)[[currSubject]])))
        }
      )
    }
    tmpExpr$cellSubset = temp_cellSubsets
    mergedExpr = rbind(mergedExpr, tmpExpr)
  }

  # building a discovr experiment S3 object
  exptInProgress <- structure(list(), class = "discovrExperiment")
  exptInProgress$markerInfoFile     <- markerInfoFile
  exptInProgress$markerInfo         <- markerInfo
  exptInProgress$fcsInfoFile        <- fcsInfoFile
  exptInProgress$fcsInfo            <- fcsInfo
  exptInProgress$parentPopulation   <- parentPopulation
  exptInProgress$allDataTransformed <- allDataTransformed
  exptInProgress$mergedExpr         <- mergedExpr
  exptInProgress$clusteringMarkers  <- clusteringMarkers
  exptInProgress$status             <- "initialized"

  return(exptInProgress)
}

#' Perform phenograph clustering for a DISCOV-R experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment()}
#' @param method A character string indicating the clustering method to use.
#' Currently only 'phenograph' is supported as a clustering method.
#' @param seed A number to be used as the seed for pseudorandom number generation (default: 12345)
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{setupDiscovrExperiment}} \code{\link{discovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom flowCore exprs
#' @importFrom igraph graph.data.frame cluster_louvain membership modularity
#' @export
clusterDiscovrExperiment <- function(
  experiment,
  method = "phenograph",
  seed = 12345,
  verbose = TRUE
){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object.",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }

  set.seed(seed)

  #########################################################################
  # Section 2.d.i from original SOP - Phenograph clustering - takes a while
  #########################################################################
  # Clustering markers object is defined in the first set-up chunk.  Tweak which markers are included up there.

  # Set up Phenograph function to use kd tree
  findNeighbors <- function(data, k){
    nearest <- RANN::nn2(data, data, k, treetype = "kd", searchtype = "standard")
    return(nearest[[1]])
  }

  runRpheno <- function(data, k=30){
    if(is.data.frame(data))
      data <- as.matrix(data)

    if(!is.matrix(data))
      stop("Wrong input data, should be a data frame or matrix!")

    if(k<1){
      stop("k must be a positive integer!")
    }else if (k > nrow(data)-2){
      stop("k must be smaller than the total number of points!")
    }

    if(verbose){
      message("Run Rphenograph starts:","\n",
              "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
              "  -k is set to ", k)
      message("  Finding nearest neighbors...")
    }

    t1 <- system.time(neighborMatrix <- findNeighbors(data, k=k+1)[,-1])

    if(verbose){
      message(
        "DONE ~",t1[3],"s\n",
        "Compute jaccard coefficient between nearest-neighbor sets..."
      )
    }
    t2 <- system.time(links <- jaccard_coeff(neighborMatrix))

    if(verbose){
      message(
        "DONE ~",t2[3],"s\n",
        "Build undirected graph from the weighted links..."
      )
    }

    links <- links[links[,1]>0, ]
    relations <- as.data.frame(links)
    colnames(relations)<- c("from","to","weight")
    t3 <- system.time(g <- igraph::graph.data.frame(relations, directed=FALSE))

    if(verbose){
      message(
        "DONE ~",t3[3],"s\n",
        "Run louvain clustering on the graph ..."
      )
    }

    t4 <- system.time(community <- igraph::cluster_louvain(g))
    if(verbose){
      message("DONE ~",t4[3],"s\n")
    }

    if(verbose){
      message("Run Rphenograph DONE, took a total of ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
      message("  Return a community class\n  -Modularity value:", igraph::modularity(community),"\n")
      message("  -Number of clusters:", length(unique(igraph::membership(community))))

    }
    return(community)
  }

  # Run phenograph (using kd treetype) on each subject.
  # Note that here 'fcs' is the flowFrame object for a subject post-transformation
  phenographClust = function(fcs, clusteringMarkers) {
    exprsMat = as.matrix(as.data.frame(exprs(fcs))[,clusteringMarkers])
    rPhenoVect = as.numeric(igraph::membership(runRpheno(data = exprsMat)))
    return(rPhenoVect)
  }

  # Cluster data in experiment
  if (method == "phenograph") {
    experiment$mergedExpr$RPclust = unlist(
      lapply(
        experiment$allDataTransformed,
        phenographClust,
        experiment$clusteringMarkers
      )
    )
  } else {
    stop("Clustering method '", method, "' is not currently supported. Stopping...")
  }

  ###################################################################
  # Sections 2.e.i from original SOP - Summarize and save outputs
  ###################################################################
  # Calculate mean expression value of each marker for each phenograph cluster in each subject
  clusterMeans <- experiment$mergedExpr %>%
    dplyr::select(-.data$cellSubset) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(.data$samp, .data$RPclust) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::mutate(RPclust = as.character(.data$RPclust))

  ##################################
  # Calculate total mean expression for each subject
  parentMeans = experiment$mergedExpr %>%
    dplyr::select(-.data$cellSubset, -.data$RPclust) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(.data$samp) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::mutate(RPclust = "Total_Parent")

  experiment$clusterMeans <- bind_rows(clusterMeans, parentMeans)

  # Count cells of each subpopulation (eg Tmr) in each phenograph cluster (from each sample)
  clusterRarePopCts <-
    experiment$mergedExpr %>%
    dplyr::select(.data$samp, .data$cellSubset, .data$RPclust) %>%
    dplyr::group_by(.data$samp, .data$RPclust) %>%
    dplyr::summarise(Total=n(), .groups = "drop_last") %>%
    as.data.frame()

  uniqueSubsets <- unique(experiment$mergedExpr$cellSubset)

  for(currCellSubset in uniqueSubsets){
    if(verbose){message("Counting cluster events for ", currCellSubset)}
    additionalMatrix = experiment$mergedExpr %>%
      dplyr::select(.data$samp, .data$cellSubset, .data$RPclust) %>%
      group_by(.data$samp, .data$RPclust) %>%
      summarise(number = sum(.data$cellSubset == currCellSubset)) %>%
      as.data.frame()
    clusterRarePopCts <- cbind(clusterRarePopCts, additionalMatrix$number)
  }

  for(i in 1:(ncol(clusterRarePopCts)-3)){
    colnames(clusterRarePopCts)[i+3] = (uniqueSubsets)[i]
  }

  aggregateCounts = clusterRarePopCts %>%
    dplyr::select(-.data$RPclust) %>%
    group_by(.data$samp) %>%
    summarise_all(sum) %>%
    rename_at(vars(-.data$samp),function(name) paste0(name,"_tot"))

  clusterRarePopCts = left_join(clusterRarePopCts, aggregateCounts)

  for(i in 1:(length(uniqueSubsets)+1)){
    clusterRarePopCts <- cbind(
      clusterRarePopCts,
      clusterRarePopCts[,i+2]/clusterRarePopCts[,i+2+length(uniqueSubsets)+1]*100
    )
  }

  for(i in 1:(length(uniqueSubsets))){
    colnames(clusterRarePopCts)[i+3+(length(uniqueSubsets)+1)*2]= (paste0("pct_", (uniqueSubsets)[i], "_in_clust"))
  }
  colnames(clusterRarePopCts)[3+(length(uniqueSubsets)+1)*2]= "pct_Total_in_clust"

  # update experiment data
  experiment$status             <- "clustered"
  experiment$clusterMethod      <- method
  experiment$clusterRarePopCts  <- clusterRarePopCts

  return(experiment)
}
