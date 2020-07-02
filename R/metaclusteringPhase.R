## Copyright (C) 2020  Mario Rosasco and Benaroya Research Institute
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
#' markers. By default all markers will be included. (default: NA)
#' 'markerInfoFile'
#' @param pctInClusterThreshold A numeric indicating a percentage. Cluster with event
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
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import dplyr
#' @importFrom flowCore exprs
#' @importFrom igraph graph.data.frame cluster_louvain membership modularity
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
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function ",
      "and perform the initial clustering using the 'clusterDiscovrExperiment' function. "
    )
  }
  if(experiment$status != "clustered"){
    stop(
      "The experiment must have status 'clustered' in order to be ready for metaclustering. ",
      "The current experiment has a status of ", experiment$status, ". ",
      "Please make sure to create your experiment using the 'setupDiscovrExperiment' function ",
      "and perform the initial clustering using the 'clusterDiscovrExperiment' function. "
    )
  }


  #########################################################################
  # Section 3.b.i-ix from original SOP - Set up metaclustering
  #########################################################################

  # Munge data, filter out low abundance clusters
  experiment$clusterRarePopCts <-
    # clean up data
    experiment$clusterRarePopCts %>%
    dplyr::mutate(RPclust = as.character(RPclust)) %>%
    left_join(experiment$clusterMeans, by = c("samp", "RPclust")) %>%
    dplyr::filter(!samp %in% dropSamples) %>%
    # Identify low abundance clusters (containing <threshold% of events per subject)
    dplyr::group_by(samp) %>%
    dplyr::mutate(totalEvents = sum(Total)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pctCellsInClust = (Total/totalEvents)*100) %>%
    dplyr::filter(pctCellsInClust >= pctInClusterThreshold)

  # remove samples and low-abundance clusters from mean data
  experiment$clusterMeans <-
    experiment$clusterMeans %>%
    dplyr::filter(!samp %in% dropSamples) %>%
    dplyr::mutate(RPclust = as.character(RPclust)) %>%
    dplyr::filter(RPclust %in% c(experiment$clusterRarePopCts$RPclust, "Total_Parent"))

  # remove samples and low-abundance clusters from expr data
  experiment$mergedExpr <-
    experiment$mergedExpr %>%
    dplyr::filter(!samp %in% dropSamples) %>%
    dplyr::mutate(RPclust = as.character(RPclust)) %>%
    dplyr::filter(RPclust %in% experiment$clusterRarePopCts$RPclust)

  # Compute the per-subject/per-cluster standard deviation to later compute z-score
  markerStdDev <-
    experiment$mergedExpr %>%
    dplyr::select(-cellSubset) %>%
    dplyr::group_by(samp, RPclust) %>%
    dplyr::summarise_all(sd)

  # Calculate parent pop standard deviations for each sample
  parentStdDev <- experiment$mergedExpr %>%
    dplyr::select(-cellSubset, -RPclust) %>%
    dplyr::group_by(samp) %>%
    dplyr::summarise_all(sd) %>%
    dplyr::mutate(RPclust = "Total_Parent")

  # Merge per-cluster and total parent standard deviations
  markerStdDev = dplyr::bind_rows(markerStdDev, parentStdDev)

  # format mean/var data for ease of use
  subjectMeans = experiment$clusterMeans %>%
    dplyr::filter(RPclust == "Total_Parent") %>%
    dplyr::select(-RPclust) %>%
    tidyr::gather("marker", "subjectMean", -samp) %>%
    dplyr::rename(subject = samp)

  subjectStdDevs = markerStdDev %>%
    dplyr::filter(RPclust == "Total_Parent") %>%
    dplyr::select(-RPclust) %>%
    gather("marker", "subjectStdDev", -samp) %>%
    dplyr::rename(subject = samp)

  subjectMeanVar <-
    dplyr::left_join(subjectMeans, subjectStdDevs, by = c("subject", "marker")) %>%
    dplyr::mutate(subjectVar = subjectStdDev **2)

  #########################################################################
  # Section 3.c from original SOP - calculate z-scores
  #########################################################################
  ##### lines 257-299 - not mentioned in SOP, but sets up data for next section

  # drop any markers that won't be used in metaclustering or plotting
  markers <- experiment$markerInfo$commonMarkerName
  markers <- markers[!markers %in% dropMarkers]

  # prepare data for each subset (ie: tmr)
  subsets <- unique(experiment$fcsInfo$cellSubset)
  dfAllSubsets <- NULL

  # build dataset of unique clusters with event occupancy for each subset
  for (currSubset in subsets){
    dfAllSubsets <- dplyr::bind_rows(
      dfAllSubsets,
      experiment$clusterRarePopCts %>%
        dplyr::rename(subject = samp) %>%
        dplyr::filter(!!as.symbol(currSubset) > 0)
    ) %>%
    unique()
  }


  # MGR - not sure what this comment means? return to this at some pt
  ## Include regardless of colorbar visualization
  clustSigPass <-
    dfAllSubsets %>%
    dplyr::select(subject, RPclust, !!subsets) %>%
    dplyr::mutate(sample = paste0(.$subject, "_", .$RPclust))

  subsetEventCounting <- clustSigPass %>% dplyr::select(sample, !!subsets)

  # for each subset (ie: tmr) compute the fraction of all events that fall in each cluster
  for(currSubset in subsets){
    currTot <- paste0("total", currSubset)
    clustSigPass <-
      dplyr::group_by(clustSigPass, subject) %>%
      dplyr::mutate(!!currTot := sum(.data[[currSubset]])) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!currSubset := .data[[currSubset]]/.data[[currTot]])
  }

  # round the fractional occupancy to 2 sig. digits and extract data w/ labels
  clustSigPass[,subsets] <- sapply(clustSigPass[,subsets], round, digits = 2)
  clustSigPass <- clustSigPass[,c("sample", subsets)]

  # Create main data object for plotting heatmaps
  hmapDfAllSubsets <-
    dfAllSubsets %>%
    dplyr::select(subject, RPclust, one_of(markers)) %>%
    reshape2::melt(id.vars = c("subject", "RPclust")) %>%
    dplyr::rename(marker = variable, mean = value) %>%
    merge(subjectMeanVar, by=c("subject", "marker")) %>%
    ## Apply z-score here!
    mutate(
      subjectZScore = (mean - subjectMean)/subjectStdDev,
      subjectArcsinhMean = mean,
      sample = paste0(.$subject, "_", .$RPclust)
    ) %>%
    dplyr::select(sample, subject, marker, subjectZScore, subjectArcsinhMean) %>%
    merge(clustSigPass, by="sample")

  # Extract the z-scores for all subjects
  allSubsetAllSubjectZscores <-
    hmapDfAllSubsets %>%
    reshape2::dcast(marker ~ sample, value.var = "subjectZScore") %>%
    magrittr::set_rownames(.$marker) %>%
    dplyr::select(-marker)

  # Extract the arcsinh fluorescence values for all subjects
  allSubsetAllSubjectArcsinh <-
    hmapDfAllSubsets %>%
    reshape2::dcast(marker ~ sample, value.var = "subjectArcsinhMean") %>%
    magrittr::set_rownames(.$marker) %>%
    dplyr::select(-marker)

  # use ComplexHeatmap as a convenient way of applying metaclustering
  metaxHeatmap <- ComplexHeatmap::Heatmap(
    allSubsetAllSubjectZscores,
    clustering_method_columns = linkage,
    clustering_distance_rows = distance,
    clustering_method_rows = linkage,
    clustering_distance_rows = distance
  )

  # cut the heatmap dendrogram to get phenotypic metaclusters
  colIndices <- stats::cutree(
    as.hclust(ComplexHeatmap::column_dend(metaxHeatmap)),
    k = nMetaclusters
  )

  ####################################################
  # attach data and return metaclustered experiment
  ####################################################
  experiment$subjectMeanVar               <- subjectMeanVar
  experiment$dfAllSubsets                 <- dfAllSubsets
  experiment$clustSigPass                 <- clustSigPass
  experiment$hmapDfAllSubsets             <- hmapDfAllSubsets
  experiment$allSubsetAllSubjectZscores   <- allSubsetAllSubjectZscores
  experiment$allSubsetAllSubjectArcsinh   <- allSubsetAllSubjectArcsinh
  experiment$subsetEventCounting          <- subsetEventCounting
  experiment$nMetaclusters                <- nMetaclusters
  experiment$colIndices                   <- colIndices
  experiment$linkage                      <- linkage
  experiment$distance                     <- distance
  experiment$status                       <- "metaclustered"

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
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import dplyr
#' @importFrom flowCore exprs
#' @importFrom igraph graph.data.frame cluster_louvain membership modularity
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
      "The experiment must have already been metaclustered in order to be re-cut using this function.. ",
      "The current experiment has a status of ", experiment$status, ". ",
      "Please use the function 'metaclusterDiscovrExperiment' if performing metaclustering for the first time."
    )
  }

  # use ComplexHeatmap as a convenient way of applying metaclustering
  metaxHeatmap <- ComplexHeatmap::Heatmap(
    allSubsetAllSubjectZscores,
    clustering_method_columns = linkage,
    clustering_distance_rows = distance,
    clustering_method_rows = linkage,
    clustering_distance_rows = distance
  )

  # cut the heatmap dendrogram to get phenotypic metaclusters
  colIndices <- stats::cutree(
    as.hclust(ComplexHeatmap::column_dend(metaxHeatmap)),
    k = nMetaclusters
  )

  experiment$nMetaclusters                <- nMetaclusters
  experiment$colIndices                   <- colIndices
  experiment$linkage                      <- linkage
  experiment$distance                     <- distance
  return(experiment)
}
