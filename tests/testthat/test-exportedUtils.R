# Input-validation tests for exported utilities (R/exportedUtils.R), using small
# mock objects. Adapted from the pilot test suite in BenaroyaResearch/briDiscovr-2,
# with error messages updated to match this version of the package.
# Behavior on real data is tested in test-eventAccessors.R, test-fcsUtilities.R,
# test-normalizeDiscovrExperiment.R, and test-umapAndDensityPlots.R.

# getSubjectCounts ------------------------------------------------------------

test_that("getSubjectCounts returns the number of events per sample", {
  result <- getSubjectCounts(create_mock_experiment())
  expect_s3_class(result, "data.frame")
  expect_named(result, c("sample", "nEvents"))
  expect_identical(result$sample, c("Sample1", "Sample2"))
  expect_equal(result$nEvents, c(50, 50))
})

test_that("getSubjectCounts validates its input", {
  expect_error(getSubjectCounts(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(
    getSubjectCounts(create_mock_experiment(include_mergedExpr = FALSE)),
    "appears incorrectly formatted"
  )
  experiment <- create_mock_experiment()
  experiment$mergedExpr$samp <- NULL
  expect_error(getSubjectCounts(experiment), "appears incorrectly formatted")
})

# getSubjectClusters ----------------------------------------------------------

test_that("getSubjectClusters returns the number of clusters per sample", {
  result <- getSubjectClusters(create_mock_experiment(status = "clustered"))
  expect_named(result, c("sample", "kClusters"))
  # the mock assigns clusters 1-3 to Sample1 and clusters 3-5 to Sample2
  expect_equal(result$kClusters, c(3, 5))
})

test_that("getSubjectClusters validates its input", {
  expect_error(getSubjectClusters(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(
    getSubjectClusters(create_mock_experiment(include_mergedExpr = FALSE)),
    "appears incorrectly formatted"
  )
  experiment <- create_mock_experiment()
  experiment$mergedExpr$RPclust <- NULL
  expect_error(getSubjectClusters(experiment), "does not appear to have been clustered yet")
})

# getMetaclusterOccupancy -----------------------------------------------------

mockOccupancyExperiment <- function(){
  experiment <- create_mock_experiment(status = "metaclustered")
  experiment$metaclusterOccupancy <- data.frame(
    subject = c("S1", "S2"),
    metacluster = c("metacluster_1", "metacluster_2"),
    Parent_frac = c(0.12345, 0.67891),
    Subset_frac = c(0.5555, 0.4444),
    Parent_cts = c(10, 20)
  )
  experiment
}

test_that("getMetaclusterOccupancy rounds only the fraction columns", {
  result <- getMetaclusterOccupancy(mockOccupancyExperiment(), precision = 2)
  expect_equal(result$Parent_frac, c(0.12, 0.68))
  expect_equal(result$Subset_frac, c(0.56, 0.44))
  expect_equal(result$Parent_cts, c(10, 20))
  expect_equal(getMetaclusterOccupancy(mockOccupancyExperiment(), precision = 3)$Parent_frac, c(0.123, 0.679))
})

test_that("getMetaclusterOccupancy validates its input", {
  experiment <- mockOccupancyExperiment()
  expect_error(getMetaclusterOccupancy(experiment, precision = 2.5), "'precision' argument must be set to an integer")
  expect_error(getMetaclusterOccupancy(experiment, precision = "two"), "'precision' argument must be set to an integer")
  expect_error(getMetaclusterOccupancy(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(getMetaclusterOccupancy(create_mock_experiment()), "does not appear to have been metaclustered yet")
})

# getEventIntensities and getEventZScores ------------------------------------

for(accessor in c("getEventIntensities", "getEventZScores")){
  accessorFunction <- get(accessor)

  test_that(paste(accessor, "validates the experiment"), {
    expect_error(accessorFunction(create_non_experiment()), "not a valid DISCOV-R experiment object")
    expect_error(
      accessorFunction(create_mock_experiment(include_mergedExpr = FALSE)),
      "appears incorrectly formatted"
    )
    unclustered <- create_mock_experiment()
    unclustered$mergedExpr$RPclust <- NULL
    expect_error(accessorFunction(unclustered), "does not appear to have been clustered yet")
    expect_error(
      accessorFunction(create_mock_experiment(status = "normalized")),
      "must be metaclustered before using this function"
    )
  })

  test_that(paste(accessor, "validates its filters"), {
    experiment <- create_mock_metaclustered_experiment()
    expect_error(accessorFunction(experiment, cellSubsets = "nonexistent"), "'cellSubset' filter could not be found")
    expect_error(accessorFunction(experiment, subjects = "nonexistent"), "'subject' filter could not be found")
    expect_error(accessorFunction(experiment, markers = "nonexistent"), "'markers' filter could not be found")
    expect_error(accessorFunction(experiment, metaclusters = 1.5), "do not appear to be whole number values")
    expect_error(accessorFunction(experiment, metaclusters = 10), "greater than the number of metaclusters")
  })
}

# checkFcsFiles (behavior on real files is in test-fcsUtilities.R) -----------

test_that("checkFcsFiles requires subject, cellSubset, and filename columns", {
  expect_error(checkFcsFiles(data.frame(foo = "bar")), "must contain columns with names")
})

test_that("checkFcsFiles reports files that do not exist", {
  fcsInfo <- data.frame(subject = "S1", cellSubset = "Parent", filename = "/nonexistent/file.fcs")
  expect_error(checkFcsFiles(fcsInfo, checkMemory = FALSE), "could not be found")
})

# testNormalizationMethodByMarker --------------------------------------------

test_that("testNormalizationMethodByMarker validates its input", {
  experiment <- create_mock_experiment(status = "clustered")
  expect_error(
    testNormalizationMethodByMarker(create_non_experiment(), markers = "CD4"),
    "must be a discovrExperiment object"
  )
  expect_error(testNormalizationMethodByMarker(experiment, markers = 123), "must be a character vector")
  expect_error(
    testNormalizationMethodByMarker(experiment, markers = "NonexistentMarker"),
    "not found in the provided experiment object"
  )
  expect_error(
    testNormalizationMethodByMarker(
      experiment, markers = "CD4", normalizationMethod = c("zScore", "none"), seed = c(123, 456)
    ),
    "Multiple random number seeds can only run with a single warpSet"
  )
  expect_error(
    testNormalizationMethodByMarker(experiment, markers = "CD4", normalizationMethod = "zScore", seed = c(123, 456)),
    "Multiple random number seeds can only run with a single warpSet"
  )
})

test_that("testNormalizationMethodByMarker returns one experiment per method, for the requested markers only", {
  experiment <- create_mock_experiment(status = "clustered")
  single <- suppressMessages(testNormalizationMethodByMarker(experiment, markers = "CD4", normalizationMethod = "zScore"))
  expect_true(is.discovrExperiment(single))
  expect_identical(single$markersNormalized, "CD4")
  expect_false(any(c("CD8", "CD3", "CD19") %in% colnames(single$mergedExprNormalizedScaled)))

  multiple <- suppressMessages(
    testNormalizationMethodByMarker(experiment, markers = c("CD4", "CD8"), normalizationMethod = c("none", "zScore"))
  )
  expect_named(multiple, c("none", "zScore"))
  expect_true(all(vapply(multiple, is.discovrExperiment, logical(1))))
})

# runUmapDiscovrExperiment (behavior on real data is in test-umapAndDensityPlots.R)

test_that("runUmapDiscovrExperiment validates its input", {
  experiment <- create_mock_normalized_experiment()
  expect_error(runUmapDiscovrExperiment(create_non_experiment()), "must be a discovrExperiment object")
  expect_error(
    runUmapDiscovrExperiment(experiment, downsampleFreq = -1),
    "'downsampleFreq' must consist of positive integers"
  )
  expect_error(
    runUmapDiscovrExperiment(experiment, downsampleFreq = 2.5),
    "'downsampleFreq' must consist of positive integers"
  )
  expect_error(
    runUmapDiscovrExperiment(experiment, downsampleBy = "number", downsampleNumber = -1),
    "'downsampleNumber' must consist of positive integers"
  )
  expect_error(
    runUmapDiscovrExperiment(experiment, downsampleFreq = c(wrongName = 1, otherName = 2)),
    "the names must match the values of cellSubset"
  )
  expect_error(runUmapDiscovrExperiment(experiment, umapMarkers = "NonexistentMarker"), "markers that are not present")
  expect_error(runUmapDiscovrExperiment(experiment, umapMarkers = "CD4"), "must contain at least two markers")
  expect_error(runUmapDiscovrExperiment(experiment, downsampleBy = "random"))
})

test_that("runUmapDiscovrExperiment requires normalization information for an un-normalized experiment", {
  experiment <- create_mock_experiment(status = "clustered")
  experiment$markerInfo$normalizationMethod <- NULL
  expect_error(runUmapDiscovrExperiment(experiment), "must include normalization specifications")
})

# plotDensityNormalizedExprsDiscovrExperiment ---------------------------------

test_that("plotDensityNormalizedExprsDiscovrExperiment validates its input", {
  expect_error(
    plotDensityNormalizedExprsDiscovrExperiment(create_non_experiment()),
    "must be a discovrExperiment object"
  )
  expect_error(
    plotDensityNormalizedExprsDiscovrExperiment(create_mock_experiment(status = "initialized")),
    "must have status"
  )
  experiment <- create_mock_experiment(status = "normalized")
  expect_error(
    plotDensityNormalizedExprsDiscovrExperiment(experiment),
    "does not include normalized expression values"
  )
  expect_error(
    plotDensityNormalizedExprsDiscovrExperiment(create_mock_normalized_experiment(), filenameOut = tempfile(fileext = ".txt")),
    "must have extension"
  )
})
