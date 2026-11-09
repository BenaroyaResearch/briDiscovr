# Tests for the per-subject and per-event accessors in R/exportedUtils.R
# (getSubjectCounts, getSubjectClusters, getEventIntensities, getEventZScores),
# using the example data. Input validation is in test-exportedUtils.R, and
# multi-value filters in test-multiValueFilters.R.

test_that("getSubjectCounts gives the number of events per subject in the FCS files", {
  skipIfNoExampleData()
  counts <- getSubjectCounts(getExampleExperiment("initialized"))
  fcsInfo <- exampleDataPaths()$fcsInfo
  expected <- tapply(vapply(fcsInfo$filename, getFcsNEvents, integer(1)), fcsInfo$subject, sum)
  expect_setequal(counts$sample, exampleSubjects)
  expect_equal(counts$nEvents, as.vector(expected[counts$sample]))
})

test_that("getSubjectClusters gives the number of clusters per subject, numbered consecutively", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  nClusters <- getSubjectClusters(clustered)
  expect_setequal(nClusters$sample, exampleSubjects)
  for(subject in nClusters$sample){
    labels <- clustered$mergedExpr$RPclust[clustered$mergedExpr$samp == subject]
    expect_identical(sort(unique(labels)), as.numeric(seq_len(nClusters$kClusters[nClusters$sample == subject])))
  }
  expect_false(any(is.na(clustered$mergedExpr$RPclust)))
})

# markers in markerInfo that are present in the expression data. The example
# markerInfo also lists FSC-A and FSC-H, which have no description in the FCS
# files and so are dropped by setupDiscovrExperiment().
markersInData <- function(experiment){
  intersect(experiment$markerInfo$commonMarkerName, colnames(experiment$mergedExpr))
}

test_that("getEventIntensities and getEventZScores work with default markers", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  expect_no_error(suppressMessages(getEventIntensities(metaclustered)))
  expect_no_error(suppressMessages(getEventZScores(metaclustered, subjects = "donor1")))
})

test_that("getEventIntensities returns the arcsinh values of every retained event, with metaclusters", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  markers <- markersInData(metaclustered)
  intensities <- suppressMessages(getEventIntensities(metaclustered, markers = markers))
  expect_equal(nrow(intensities), nrow(metaclustered$mergedExpr))
  expect_setequal(colnames(intensities), c(markers, "subject", "cellSubset", "cluster", "metacluster"))
  expect_equal(intensities$CD20, metaclustered$mergedExpr$CD20)
  expect_equal(intensities$metacluster, unname(metaclustered$colIndices[intensities$cluster]))
})

test_that("getEventIntensities filters return only the requested events and markers", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  result <- suppressMessages(getEventIntensities(
    metaclustered, cellSubsets = exampleRareSubset, subjects = "donor4", metaclusters = 2, markers = c("CD27", "IgD")
  ))
  expect_true(all(result$subject == "donor4"))
  expect_true(all(unlist(result$cellSubset) == exampleRareSubset))
  expect_true(all(result$metacluster == 2))
  expect_setequal(colnames(result), c("CD27", "IgD", "subject", "cellSubset", "cluster", "metacluster"))
  expected <- metaclustered$mergedExpr$samp == "donor4" &
    unlist(metaclustered$mergedExpr$cellSubset) == exampleRareSubset &
    metaclustered$colIndices[metaclustered$mergedExpr$sampRpClust] == 2
  expect_equal(nrow(result), sum(expected))
  expect_equal(result$IgD, metaclustered$mergedExpr$IgD[expected])
})

test_that("getEventZScores standardizes each marker within each subject", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  zScores <- suppressMessages(getEventZScores(metaclustered, subjects = c("donor1", "donor5"), markers = markersInData(metaclustered)))
  expect_equal(nrow(zScores), sum(metaclustered$mergedExpr$samp %in% c("donor1", "donor5")))
  for(subject in c("donor1", "donor5")){
    isSubject <- metaclustered$mergedExpr$samp == subject
    for(marker in c("CD20", "CD138")){
      values <- metaclustered$mergedExpr[[marker]][isSubject]
      expected <- (values - mean(values)) / sd(values)
      observed <- zScores[[marker]][zScores$subject == subject]
      # getEventZScores reorders events, so compare the sorted values
      expect_equal(sort(observed), sort(expected))
    }
  }
})

test_that("getEventZScores uses the whole subject for the mean, even when filtering", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  rareOnly <- suppressMessages(getEventZScores(metaclustered, subjects = "donor2", cellSubsets = exampleRareSubset, markers = "CD38"))
  isSubject <- metaclustered$mergedExpr$samp == "donor2"
  isRare <- unlist(metaclustered$mergedExpr$cellSubset) == exampleRareSubset
  values <- metaclustered$mergedExpr$CD38
  expected <- (values[isSubject & isRare] - mean(values[isSubject])) / sd(values[isSubject])
  expect_equal(sort(rareOnly$CD38), sort(expected))
})
