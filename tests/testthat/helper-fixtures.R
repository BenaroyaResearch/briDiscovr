# Small hand-built objects for fast input-validation tests, without reading FCS files.
# Adapted from the pilot test suite in BenaroyaResearch/briDiscovr-2.
# The structure follows what setupDiscovrExperiment() and clusterDiscovrExperiment()
# produce: markerInfo has columns commonMarkerName/fcsMarkerName/useToCluster/
# normalizationMethod, and mergedExpr$cellSubset is a list column.

# A mock discovrExperiment with two samples ("Sample1", "Sample2") of 50 events
# each, four markers (CD4, CD8, CD3 used for clustering; CD19 not), and five
# clusters (RPclust 1-5) if include_mergedExpr is TRUE.
create_mock_experiment <- function(status = "initialized",
                                   include_mergedExpr = TRUE,
                                   include_clusterMeans = FALSE){
  markerInfo <- data.frame(
    commonMarkerName = c("CD4", "CD8", "CD3", "CD19"),
    fcsMarkerName = c("CD4_BV421", "CD8_FITC", "CD3_PE", "CD19_APC"),
    useToCluster = c(TRUE, TRUE, TRUE, FALSE),
    normalizationMethod = c("zScore", "zScore", "zScore", "none"),
    stringsAsFactors = FALSE
  )

  experiment <- list(
    status = status,
    markerInfo = markerInfo,
    parentPopulation = "Parent",
    clusteringMarkers = c("CD4", "CD8", "CD3")
  )

  if(include_mergedExpr){
    withr::with_seed(42, {
      n <- 100
      experiment$mergedExpr <- data.frame(
        CD19 = rnorm(n, 2, 1),
        CD3 = rnorm(n, 7, 1),
        CD4 = rnorm(n, 5, 1),
        CD8 = rnorm(n, 3, 1),
        samp = rep(c("Sample1", "Sample2"), each = n / 2),
        RPclust = rep(1:5, each = n / 5),
        stringsAsFactors = FALSE
      )
      experiment$mergedExpr$cellSubset <- as.list(rep(c("Parent", "Subset"), n / 2))
    })
  }

  if(include_clusterMeans){
    withr::with_seed(43, {
      experiment$clusterMeans <- data.frame(
        samp = rep(c("Sample1", "Sample2"), each = 5),
        RPclust = as.character(rep(1:5, 2)),
        CD4 = rnorm(10, 5, 0.5),
        CD8 = rnorm(10, 3, 0.5),
        CD3 = rnorm(10, 7, 0.5),
        CD19 = rnorm(10, 2, 0.5),
        stringsAsFactors = FALSE
      )
    })
  }

  class(experiment) <- "discovrExperiment"
  experiment
}

# A mock experiment that looks normalized (normalized expression = raw expression),
# enough to reach the argument checks in runUmapDiscovrExperiment().
create_mock_normalized_experiment <- function(){
  experiment <- create_mock_experiment(status = "normalized")
  experiment$mergedExprNormalizedScaled <- experiment$mergedExpr
  experiment
}

# A mock experiment that looks metaclustered, enough to reach the filter checks
# in getEventIntensities() and getEventZScores().
create_mock_metaclustered_experiment <- function(){
  experiment <- create_mock_experiment(status = "metaclustered")
  experiment$nMetaclusters <- 3
  experiment$mergedExpr$RPclust <- as.character(experiment$mergedExpr$RPclust)
  experiment$mergedExpr$sampRpClust <- paste0(experiment$mergedExpr$samp, "_", experiment$mergedExpr$RPclust)
  sampRpClusts <- unique(experiment$mergedExpr$sampRpClust)
  experiment$colIndices <- setNames(rep_len(1:3, length(sampRpClusts)), sampRpClusts)
  experiment
}

# A list that is not a discovrExperiment, for testing input validation.
create_non_experiment <- function(){
  list(foo = "bar", data = 1:10)
}
