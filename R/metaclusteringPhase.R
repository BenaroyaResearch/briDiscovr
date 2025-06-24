## Copyright (C) 2025  Matt Dufort, Mario Rosasco, Virginia Muir, and Benaroya Research Institute
##
## This file is part of the briDiscovr package

#' Perform metaclustering for a DISCOV-R experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment}
#' and clustered using \code{clusterDiscovrExperiment}
#' @param dropSamples A string or vector of strings indicating names of samples
#' to exclude from further analysis. Note that the convenience function
#' \code{getSampleNames} can be used to view a list of all sample names. (default: NA)
#' @param dropMarkers A string or vector of strings indicating names of markers to exclude
#' from analysis. Note that printing the experiment object will display a list of all
#' markers. By default all clustering markers will be included. (default: NA)
#' @param pctInClusterThreshold A numeric indicating a percentage. Clusters with event
#' occupancy below this threshold will be considered 'low-abundance' and will be
#' dropped from further analysis (default: 1)
#' @param nMetaclusters A numeric indicating the number of metaclusters to generate (default: 12)
#' @param linkage A string indicating the cluster linkage method to use when metaclustering.
#' See \code{hclust} for more details. (default: "ward.D2")
#' @param distance A string indicating the distance metric to use when metaclustering. (default: "euclidean")
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{setupDiscovrExperiment}} \code{\link{clusterDiscovrExperiment}}
#' @author Mario G Rosasco, Virginia Muir
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @import dplyr
#' @importFrom rlang sym :=
#' @importFrom tidyr gather
#' @importFrom tibble column_to_rownames
#' @importFrom flowCore exprs
#' @importFrom stats as.hclust sd cutree
#' @export
metaclusterDiscovrExperiment <- function(
  experiment,
  dropSamples = NA,
  dropMarkers = NA,
  pctInClusterThreshold = 1,
  nMetaclusters = 12,
  linkage = "ward.D2",
  distance = "euclidean",
  verbose = TRUE
){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object.\n",
      "Please create your experiment using the 'setupDiscovrExperiment' function\n",
      "and perform the initial clustering using the 'clusterDiscovrExperiment' function. "
    )
  }
  if(!(experiment$status %in% c("clustered", "normalized"))){
    stop(
      "The experiment must have status 'clustered' or 'normalized' in order to be ready for metaclustering.\n",
      "The current experiment has a status of ", experiment$status, ".\n",
      "Please make sure to create your experiment using the 'setupDiscovrExperiment' function,\n",
      "perform the initial clustering using the 'clusterDiscovrExperiment' function,\n",
      "and optionally normalize marker expression using the 'normalizeDiscovrExperiment' function."
    )
  }
  
  if(experiment$status == "clustered"){
    message(
      "The current experiment has been clustered, but has not yet had marker expression normalized using the\n",
      "'normalizeDiscovrExperiment' function. Normalization will now be performed prior to metaclustering,\n",
      "using the normalizationMethod specifications provided in 'markerInfo', or the default normalization method\n",
      "set in 'normalizeDiscovrExperiment'.\n",
      "NOTE: This is a change in workflow from briDiscovr version 0.3, which applied z-score normalization\n",
      "as part of the metaclustering process.\n\n")
    
    experiment <- normalizeDiscovrExperiment(experiment)
  }


  #########################################################################
  # Section 3.b.i-ix from original SOP - Set up metaclustering
  #########################################################################

  # Munge data, filter out low abundance clusters
  experiment$clusterRarePopCts <-
    # clean up data
    experiment$clusterRarePopCts %>%
    dplyr::mutate(RPclust = as.character(RPclust)) %>%
    dplyr::mutate(sampRpClust = paste0(samp, "_", RPclust)) %>%
    left_join(experiment$clusterMeansNormalizedScaled, by = c("samp", "RPclust")) %>%
    dplyr::filter(!(samp %in% dropSamples)) %>%
    # Identify low abundance clusters (containing < threshold % of events per subject)
    dplyr::group_by(samp) %>%
    dplyr::mutate(totalEvents = sum(Total)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pctCellsInClust = (Total / totalEvents)*100) %>%
    dplyr::filter(pctCellsInClust >= pctInClusterThreshold)

  # get total parent cluster names for each sample to enable filtering
  sampTotalParentNames <- unique(paste0(experiment$clusterMeansNormalizedScaled$samp, "_Total_Parent"))

  # remove samples to be dropped and low-abundance clusters from mean data
  experiment$clusterMeans <-
    experiment$clusterMeans %>%
    dplyr::filter(!samp %in% dropSamples) %>%
    dplyr::mutate(
      RPclust = as.character(RPclust),
      sampRpClust = paste0(samp, "_", RPclust)
    ) %>%
    dplyr::filter(sampRpClust %in% c(experiment$clusterRarePopCts$sampRpClust, sampTotalParentNames))
  experiment$clusterMeansNormalizedScaled <-
    experiment$clusterMeansNormalizedScaled %>%
    dplyr::filter(!samp %in% dropSamples) %>%
    dplyr::mutate(
      RPclust = as.character(RPclust),
      sampRpClust = paste0(samp, "_", RPclust)
    ) %>%
    dplyr::filter(sampRpClust %in% c(experiment$clusterRarePopCts$sampRpClust, sampTotalParentNames))

  # remove samples to be dropped and low-abundance clusters from expr data
  experiment$mergedExpr <-
    experiment$mergedExpr %>%
    dplyr::filter(!(samp %in% dropSamples)) %>%
    dplyr::mutate(
      RPclust = as.character(RPclust),
      sampRpClust = paste0(samp, "_", RPclust)
    ) %>%
    dplyr::filter(sampRpClust %in% experiment$clusterRarePopCts$sampRpClust)
  experiment$mergedExprNormalizedScaled <-
    experiment$mergedExprNormalizedScaled %>%
    dplyr::filter(!(samp %in% dropSamples)) %>%
    dplyr::mutate(
      RPclust = as.character(RPclust),
      sampRpClust = paste0(samp, "_", RPclust)
    ) %>%
    dplyr::filter(sampRpClust %in% experiment$clusterRarePopCts$sampRpClust)

  # format means for ease of use
  subjectMeansNormalizedScaled <- experiment$clusterMeansNormalizedScaled %>%
    dplyr::filter(RPclust == "Total_Parent") %>%
    dplyr::select(-RPclust, -sampRpClust) %>%
    tidyr::gather("marker", "subjectMean", -samp) %>%
    dplyr::rename(subject = samp)

  #########################################################################
  # Section 3.c from original SOP - calculate z-scores
  #########################################################################
  ##### lines 257-299 - not mentioned in SOP, but sets up data for next section

  # get list of markers to be used in metaclustering
  metaclusterMarkers <- experiment$clusteringMarkers[!experiment$clusteringMarkers %in% dropMarkers]

  # prepare data for each subset (ie: tmr)
  cellSubsets <- unique(experiment$fcsInfo$cellSubset)
  dfAllSubsets <- NULL

  # build dataset of unique clusters with event occupancy for each subset
  for (currCellSubset in cellSubsets){
    dfAllSubsets <- dplyr::bind_rows(
      dfAllSubsets,
      experiment$clusterRarePopCts %>%
        dplyr::rename(subject = samp) %>%
        dplyr::filter(!!as.symbol(currCellSubset) > 0)
    ) %>%
    distinct()
  }

  # Get cluster occupancy data in easy to use format
  clustSigPass <-
    dfAllSubsets %>%
    dplyr::select(subject, RPclust, !!cellSubsets) %>%
    dplyr::mutate(sampRpClust = paste0(subject, "_", RPclust))

  subsetEventCounting <- clustSigPass %>% dplyr::select(sampRpClust, !!cellSubsets)

  # for each subset (ie: tmr) compute the fraction of all events that fall in each cluster
  for(currCellSubset in cellSubsets){
    currTot <- paste0("total", currCellSubset)
    clustSigPass <-
      clustSigPass %>%
      dplyr::group_by(subject) %>%
      dplyr::mutate(!!currTot := sum(!!rlang::sym(currCellSubset))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!currCellSubset := !!rlang::sym(currCellSubset) / !!rlang::sym(currTot))
  }

  # round the fractional occupancy to 2 sig. digits and extract data w/ labels
  clustSigPass[,cellSubsets] <- sapply(clustSigPass[,cellSubsets], round, digits = 2)
  clustSigPass <- clustSigPass[,c("sampRpClust", cellSubsets)]

  # Create main data object for plotting heatmaps
  hmapDfAllSubsets <-
    dfAllSubsets %>%
    dplyr::select(subject, RPclust, any_of(experiment$markerInfo$commonMarkerName)) %>%
    reshape2::melt(id.vars = c("subject", "RPclust")) %>%
    dplyr::rename(marker = variable, mean = value) %>%
    merge(subjectMeansNormalizedScaled, by = c("subject", "marker")) %>%
    ## Apply z-score here!
    mutate(
      subjectNormalizedScaled = mean,
      sampRpClust = paste0(subject, "_", RPclust)
    ) %>%
    dplyr::select(
      sampRpClust,
      subject,
      marker,
      subjectNormalizedScaled
    ) %>%
    merge(clustSigPass, by="sampRpClust")

  # Extract the normalized scaled means for all subjects
  allSubsetAllSubjectNormalizedScaled <-
    hmapDfAllSubsets %>%
    reshape2::dcast(marker ~ sampRpClust, value.var = "subjectNormalizedScaled") %>%
    tibble::column_to_rownames("marker")

  # use ComplexHeatmap as a convenient way of applying metaclustering
  metaxHeatmap <- 
    suppressMessages(  # suppress default messages from Heatmap
      ComplexHeatmap::Heatmap(
        as.matrix(allSubsetAllSubjectNormalizedScaled[metaclusterMarkers,]),
        clustering_method_columns = linkage,
        clustering_distance_columns = distance,
        clustering_method_rows = linkage,
        clustering_distance_rows = distance
      )
    )
  
  # cut the heatmap dendrogram to get phenotypic metaclusters
  clusterDendrogram <-
    suppressWarnings( # suppress warning messages from ComplexHeatmap
      as.hclust(ComplexHeatmap::column_dend(metaxHeatmap))
    )
  colIndices <- cutree(
    clusterDendrogram,
    k = nMetaclusters
  )
  # actual number of groups after cutting the tree
  kGroups <- length(unique(colIndices))

  # compute fraction of each cell subset in each metacluster for each sample
  metaclusterOccupancy <-
    data.frame(
      sampRpClust = names(colIndices),
      metacluster = paste0("metacluster_", colIndices)
    ) %>%
    dplyr::left_join(subsetEventCounting, by = "sampRpClust") %>%
    dplyr::mutate(subject = stringr::str_remove(sampRpClust, "_[0-9]+$"))

  # for each subset (ie: tmr) compute the fraction of all events that fall in each cluster
  for(currCellSubset in cellSubsets){
    currTot <- paste0("total", currCellSubset)
    metaclusterOccupancy <-
      dplyr::group_by(metaclusterOccupancy, subject) %>%
      dplyr::mutate(!!currTot := sum(!!rlang::sym(currCellSubset))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(subject, metacluster) %>%
      dplyr::mutate(
        !!paste0(currCellSubset, "_frac") := sum(!!rlang::sym(currCellSubset))/!!rlang::sym(currTot),
        !!paste0(currCellSubset, "_cts") := sum(!!rlang::sym(currCellSubset))
      ) %>%
      dplyr::ungroup()
  }

  # extract relevant data w/ labels
  cellsSubsetColnames <- c(paste0(cellSubsets, "_frac"), paste0(cellSubsets, "_cts"))
  metaclusterOccupancy <- metaclusterOccupancy[,c("subject", "metacluster", cellsSubsetColnames)]
  metaclusterOccupancy <- distinct(metaclusterOccupancy)

  ####################################################
  # attach data and return metaclustered experiment
  ####################################################
  experiment$status                               <- "metaclustered"
  # experiment$subjectMeanVar                       <- subjectMeanVar
  experiment$dfAllSubsets                         <- dfAllSubsets
  experiment$clustSigPass                         <- clustSigPass
  experiment$hmapDfAllSubsets                     <- hmapDfAllSubsets
  experiment$pctInClusterThreshold                <- pctInClusterThreshold
  experiment$allSubsetAllSubjectNormalizedScaled  <- allSubsetAllSubjectNormalizedScaled
  # experiment$allSubsetAllSubjectArcsinh           <- allSubsetAllSubjectArcsinh
  experiment$subsetEventCounting                  <- subsetEventCounting
  experiment$metaclusterMarkers                   <- metaclusterMarkers
  experiment$nMetaclusters                        <- nMetaclusters
  experiment$metaclusterOccupancy                 <- metaclusterOccupancy
  experiment$colIndices                           <- colIndices
  experiment$kGroups                              <- kGroups
  experiment$linkage                              <- linkage
  experiment$distance                             <- distance
  experiment$clusterDendrogram                    <- clusterDendrogram

  return(experiment)
}


#' Re-cut a metaclustered experiment with a different target number of metaclusters
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment},
#' clustered using \code{clusterDiscovrExperiment}, and metaclustered using
#' \code{metaclusterDiscovrExperiment}
#' @param nMetaclusters A numeric indicating the number of metaclusters to generate (default: 12)
#' @param linkage A string indicating the cluster linkage method to use when metaclustering.
#' See \code{hclust} for more details. (default: "ward.D2")
#' @param distance A string indicating the distance metric to use when metaclustering. (default: "euclidean")
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{metaclusterDiscovrExperiment}}
#' @author Mario G Rosasco, Virginia Muir
#' @import dplyr
#' @importFrom flowCore exprs
#' @importFrom stats as.hclust sd cutree
#' @export
recutMetaclusters <- function(
  experiment,
  nMetaclusters = 12,
  linkage = "ward.D2",
  distance = "euclidean",
  verbose = TRUE
){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function ",
      "and perform the initial clustering using the 'clusterDiscovrExperiment' function. "
    )
  }
  if(experiment$status != "metaclustered"){
    stop(
      "The experiment must have already been metaclustered in order to be re-cut using this function. ",
      "The current experiment has a status of ", experiment$status, ". ",
      "Please use the function 'metaclusterDiscovrExperiment' if performing metaclustering for the first time."
    )
  }
  if(!"clusterDendrogram" %in% names(experiment)){
    stop(
      "The object does not have the data expected in a metaclustered experiment. ",
      "Please re-run the clustering and metaclustering steps and try again."
    )
  }

  # use ComplexHeatmap as a convenient way of applying metaclustering
  metaxHeatmap <-
    suppressMessages(  # suppress default messages from Heatmap
      ComplexHeatmap::Heatmap(
        as.matrix(experiment$allSubsetAllSubjectNormalizedScaled[experiment$metaclusterMarkers,]),
        clustering_method_columns = linkage,
        clustering_distance_columns = distance,
        clustering_method_rows = linkage,
        clustering_distance_rows = distance
      )
    )
  
  # cut the heatmap dendrogram to get phenotypic metaclusters
  clusterDendrogram <-
    suppressWarnings( # suppress warning messages from ComplexHeatmap
      as.hclust(ComplexHeatmap::column_dend(metaxHeatmap))
    )
  colIndices <- cutree(
    clusterDendrogram,
    k = nMetaclusters
  )
  # actual number of groups after cutting the tree
  kGroups <- length(unique(colIndices))

  # compute fraction of each cell subset in each metacluster for each sample
  metaclusterOccupancy <-
    data.frame(
      sampRpClust = names(colIndices),
      metacluster = paste0("metacluster_", colIndices)
    ) %>%
    dplyr::left_join(experiment$subsetEventCounting, by = "sampRpClust") %>%
    dplyr::mutate(subject = stringr::str_remove(sampRpClust, "_[0-9]+$"))

  # for each subset (ie: tmr) compute the fraction of all events that fall in each cluster
  cellSubsets <- unique(experiment$fcsInfo$cellSubset)

  for(currCellSubset in cellSubsets){
    currTot <- paste0("total", currCellSubset)
    metaclusterOccupancy <-
      dplyr::group_by(metaclusterOccupancy, subject) %>%
      dplyr::mutate(!!currTot := sum(.data[[currCellSubset]])) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(subject, metacluster) %>%
      dplyr::mutate(
        !!paste0(currCellSubset, "_frac") := sum(.data[[currCellSubset]])/.data[[currTot]],
        !!paste0(currCellSubset, "_cts") := sum(.data[[currCellSubset]])
      ) %>%
      dplyr::ungroup()
  }

  # extract relevant data w/ labels
  cellsSubsetColnames = c(paste0(cellSubsets, "_frac"), paste0(cellSubsets, "_cts"))
  metaclusterOccupancy <- metaclusterOccupancy[,c("subject", "metacluster", cellsSubsetColnames)]
  metaclusterOccupancy <- unique(metaclusterOccupancy)

  experiment$nMetaclusters                <- nMetaclusters
  experiment$metaclusterOccupancy         <- metaclusterOccupancy
  experiment$colIndices                   <- colIndices
  experiment$kGroups                      <- kGroups
  experiment$linkage                      <- linkage
  experiment$distance                     <- distance
  experiment$clusterDendrogram            <- clusterDendrogram
  return(experiment)
}
