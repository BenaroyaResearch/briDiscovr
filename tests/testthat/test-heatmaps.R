# Tests for makeMetaclusterHeatmaps() (R/heatmapGenerators.R), using the example data.
# These check that the expected image files are written and that inputs are
# validated; they do not compare image contents.

expectedHeatmapFiles <- function(prefix, experiment, parentTitle = "Parent"){
  fullPrefix <- paste0(
    prefix, experiment$linkage, "_", experiment$distance, "_",
    experiment$pctInClusterThreshold, "pctThresh_", experiment$nMetaclusters, "metaclusters"
  )
  paste0(fullPrefix, c(
    "_allClusters_normCount.png",
    "_allClusters_metacluster_MFI.png",
    "_metaclusters_arcsinh_weightedAvgByClust.png",
    paste0("_total", parentTitle, "_arcsinh_weightedAvg.png"),
    # the example data has non-clustering markers (CD3, CD14), so these are written too
    "_allClusters_nonClusteringMarkers_MFI.png",
    paste0("_total", parentTitle, "_nonClusteringMarkers_arcsinh_weightedAvg.png"),
    "_nonClusteringMarkers_arcsinh_weightedAvgByClust.png"
  ))
}

test_that("makeMetaclusterHeatmaps writes the expected PNG files", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  outDir <- withr::local_tempdir()
  withr::local_dir(outDir)
  suppressMessages(makeMetaclusterHeatmaps(metaclustered, filenamePrefix = "test_", parentTitle = "TetNeg", verbose = FALSE))
  expected <- expectedHeatmapFiles("test_", metaclustered, parentTitle = "TetNeg")
  expect_setequal(list.files(outDir), expected)
  # each file is a non-empty PNG image
  for(file in expected){
    expect_identical(readBin(file, "raw", 4), as.raw(c(0x89, 0x50, 0x4e, 0x47)), info = file)
  }
})

test_that("makeMetaclusterHeatmaps uses today's date as the default file name prefix", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  outDir <- withr::local_tempdir()
  withr::local_dir(outDir)
  suppressMessages(makeMetaclusterHeatmaps(metaclustered, verbose = FALSE))
  expect_true(all(startsWith(list.files(outDir), format(Sys.Date(), "%y%m%d-"))))
})

test_that("makeMetaclusterHeatmaps accepts custom colors and dropped non-clustering markers", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  outDir <- withr::local_tempdir()
  withr::local_dir(outDir)
  expect_no_error(suppressMessages(makeMetaclusterHeatmaps(
    metaclustered, filenamePrefix = "custom_", dropMarkers = "CD14",
    metaclusterColors = grDevices::hcl.colors(exampleNMetaclusters),
    normCountBreaks = c(-3, 0, 3), intensityBreaks = c(0, 3, 6, 9),
    verbose = FALSE
  )))
  expect_length(list.files(outDir), 7)
})

test_that("makeMetaclusterHeatmaps validates its input", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  outDir <- withr::local_tempdir()
  withr::local_dir(outDir)
  expect_error(makeMetaclusterHeatmaps(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(
    makeMetaclusterHeatmaps(getExampleExperiment("normalized")),
    "must have already been metaclustered"
  )
  incomplete <- metaclustered
  incomplete$hmapDfAllSubsets <- NULL
  expect_error(makeMetaclusterHeatmaps(incomplete), "does not have the data expected")
  expect_error(
    makeMetaclusterHeatmaps(metaclustered, dropMarkers = "NotAMarker", verbose = FALSE),
    "markers were not found in the dataset: NotAMarker"
  )
  expect_error(
    makeMetaclusterHeatmaps(metaclustered, dropMarkers = "CD20", verbose = FALSE),
    "were used for metaclustering and cannot be dropped"
  )
  expect_error(
    suppressMessages(makeMetaclusterHeatmaps(metaclustered, childSubsets = "CD8", verbose = FALSE)),
    "requested child subsets.*was not found"
  )
  expect_error(
    suppressMessages(makeMetaclusterHeatmaps(metaclustered, metaclusterColors = c("red", "blue"), verbose = FALSE)),
    "length of 'metaclusterColors' must be the same as the number of metaclusters"
  )
  expect_error(
    suppressMessages(makeMetaclusterHeatmaps(metaclustered, normCountBreaks = c(-1, 1), verbose = FALSE)),
    "normCountBreaks and normCountColors must be the same length"
  )
  expect_error(
    suppressMessages(makeMetaclusterHeatmaps(metaclustered, intensityBreaks = 1:3, verbose = FALSE)),
    "intensityBreaks and intensityColors must be the same length"
  )
  # no files are written when the input is invalid
  expect_length(list.files(outDir), 0)
})
