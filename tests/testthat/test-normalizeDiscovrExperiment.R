# Tests for normalizeDiscovrExperiment() (R/clusteringPhase.R) and
# testNormalizationMethodByMarker() (R/exportedUtils.R), using the example data.
#
# Normalization has two steps: (1) each marker is normalized within each sample by
# its method (zScore, warpSet[#], or none), then (2) every normalized marker is
# scaled across all events to mean 0 and standard deviation 1.

markersWithMethod <- function(experiment, method){
  experiment$markerInfo$commonMarkerName[experiment$markerInfo$normalizationMethod == method]
}

test_that("normalization sets status and records which markers were normalized", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  expect_identical(normalized$status, "normalized")
  markerColumns <- setdiff(colnames(getExampleExperiment("clustered")$mergedExpr), c("samp", "cellSubset", "RPclust"))
  expect_setequal(normalized$markersNormalized, markerColumns)
  # markers in markerInfo that are not in the expression data (FSC) are marked as not normalized
  notInData <- setdiff(normalized$markerInfo$commonMarkerName, markerColumns)
  expect_true(all(normalized$markerInfo$normalizationMethod[normalized$markerInfo$commonMarkerName %in% notInData] == "not normalized"))
  # methods from the markerInfo file are kept for the others
  expect_setequal(markersWithMethod(normalized, "none"), c("IgM", "CD27"))
  expect_setequal(markersWithMethod(normalized, "warpSet3"), "IgD")
})

test_that("normalization keeps the arcsinh values and clustering results unchanged", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  normalized <- getExampleExperiment("normalized")
  expect_identical(normalized$mergedExpr, clustered$mergedExpr)
  expect_identical(normalized$mergedExprNormalizedScaled$RPclust, clustered$mergedExpr$RPclust)
  expect_identical(normalized$mergedExprNormalizedScaled$samp, clustered$mergedExpr$samp)
})

test_that("every normalized marker has mean 0 and standard deviation 1 across all events", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  values <- normalized$mergedExprNormalizedScaled[, normalized$markersNormalized]
  expect_equal(unname(colMeans(values)), rep(0, ncol(values)), tolerance = 1e-8)
  expect_equal(unname(apply(values, 2, sd)), rep(1, ncol(values)), tolerance = 1e-8)
})

test_that("zScore markers have mean 0 and equal spread within every sample", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  for(marker in markersWithMethod(normalized, "zScore")){
    bySample <- split(normalized$mergedExprNormalizedScaled[[marker]], normalized$mergedExprNormalizedScaled$samp)
    expect_equal(unname(vapply(bySample, mean, numeric(1))), rep(0, length(bySample)), tolerance = 1e-8)
    sds <- vapply(bySample, sd, numeric(1))
    expect_equal(unname(sds), rep(sds[[1]], length(sds)), tolerance = 1e-8)
  }
})

test_that("markers with method 'none' are only scaled across all events", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  normalized <- getExampleExperiment("normalized")
  for(marker in markersWithMethod(normalized, "none")){
    expect_equal(normalized$mergedExprNormalizedScaled[[marker]], as.vector(scale(clustered$mergedExpr[[marker]])))
  }
})

test_that("warpSet changes the values of warpSet markers", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  normalized <- getExampleExperiment("normalized")
  for(marker in c(markersWithMethod(normalized, "warpSet"), "IgD")){
    expect_false(isTRUE(all.equal(
      normalized$mergedExprNormalizedScaled[[marker]], as.vector(scale(clustered$mergedExpr[[marker]]))
    )))
    # the order of events within each sample is kept (warpSet is monotonic within a sample)
    for(sample in c("donor1", "donor8")){
      isSample <- clustered$mergedExpr$samp == sample
      expect_gt(
        cor(normalized$mergedExprNormalizedScaled[[marker]][isSample], clustered$mergedExpr[[marker]][isSample], method = "spearman"),
        0.99
      )
    }
  }
})

test_that("normalized cluster means are the means of normalized events, with a Total_Parent row per sample", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  means <- normalized$clusterMeansNormalizedScaled
  expect_setequal(means$samp[means$RPclust == "Total_Parent"], exampleSubjects)
  events <- normalized$mergedExprNormalizedScaled
  isCluster <- events$samp == "donor2" & events$RPclust == 3
  expect_equal(
    means$CD20[means$samp == "donor2" & means$RPclust == "3"],
    mean(events$CD20[isCluster])
  )
  expect_equal(
    means$CD20[means$samp == "donor2" & means$RPclust == "Total_Parent"],
    mean(events$CD20[events$samp == "donor2"])
  )
})

test_that("normalizationInfo overrides markerInfo, with the default for unlisted markers", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  normalized <- suppressMessages(normalizeDiscovrExperiment(
    clustered, normalizationInfo = c(CD20 = "none"), defaultNormalizationMethod = "zScore", verbose = FALSE
  ))
  expect_identical(normalized$markerInfo$normalizationMethod[normalized$markerInfo$commonMarkerName == "CD20"], "none")
  expect_equal(normalized$mergedExprNormalizedScaled$CD20, as.vector(scale(clustered$mergedExpr$CD20)))
  expect_true(all(markersWithMethod(normalized, "zScore") %in% setdiff(normalized$markersNormalized, "CD20")))

  allNone <- suppressMessages(normalizeDiscovrExperiment(clustered, normalizationInfo = "none", verbose = FALSE))
  expect_equal(allNone$mergedExprNormalizedScaled$IgD, as.vector(scale(clustered$mergedExpr$IgD)))
})

test_that("re-normalizing a normalized experiment gives a message and overwrites the results", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  expect_message(
    renormalized <- normalizeDiscovrExperiment(normalized, normalizationInfo = "none", verbose = FALSE),
    "normalization has already been run"
  )
  expect_equal(renormalized$mergedExprNormalizedScaled$IgD, as.vector(scale(normalized$mergedExpr$IgD)))
})

test_that("normalizeDiscovrExperiment validates its input", {
  expect_error(normalizeDiscovrExperiment(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(
    normalizeDiscovrExperiment(create_mock_experiment(status = "initialized")),
    "must have status 'clustered', 'normalized', or 'metaclustered'"
  )
  clustered <- create_mock_experiment(status = "clustered")
  expect_error(
    normalizeDiscovrExperiment(clustered, defaultNormalizationMethod = "quantile", verbose = FALSE),
    "must be 'zScore', 'none', or 'warpset\\[#\\]'"
  )
  expect_error(
    normalizeDiscovrExperiment(clustered, normalizationInfo = matrix("none", 1, 1, dimnames = list(NULL, "CD4")), verbose = FALSE),
    "Format of 'normalizationInfo' not recognized"
  )
  expect_error(
    normalizeDiscovrExperiment(clustered, normalizationInfo = c(CD4 = "quantile"), verbose = FALSE),
    "can only contain values"
  )
})

test_that("testNormalizationMethodByMarker gives reproducible warpSet results for each seed", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  results <- quietly(testNormalizationMethodByMarker(
    clustered, markers = "CD38", normalizationMethod = "warpSet", seed = c(1, 2)
  ))
  expect_named(results, c("warpSet_seed1", "warpSet_seed2"))
  again <- quietly(testNormalizationMethodByMarker(clustered, markers = "CD38", normalizationMethod = "warpSet", seed = 1))
  expect_identical(again$mergedExprNormalizedScaled$CD38, results$warpSet_seed1$mergedExprNormalizedScaled$CD38)
})
