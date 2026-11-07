# Regression tests: filter arguments that default to NA must also accept vectors
# with more than one value. Previously these used `if(is.na(x))`, which errors in
# R >= 4.2 ("the condition has length > 1") when x has length > 1.

# one small metaclustered experiment, shared by the tests in this file
metaclusteredSynthetic <- local({
  cached <- NULL
  function(){
    if(is.null(cached)){
      experiment <- makeSyntheticExperiment(nSubjects = 3, rareSubsetNames = c("rareA", "rareB"))
      experiment <- suppressMessages(clusterDiscovrExperiment(experiment, k = 10, nThreads = 1, verbose = FALSE))
      experiment <- suppressMessages(normalizeDiscovrExperiment(experiment, verbose = FALSE))
      cached <<- suppressMessages(metaclusterDiscovrExperiment(experiment, nMetaclusters = 3, verbose = FALSE))
    }
    cached
  }
})

for(accessor in c("getEventIntensities", "getEventZScores")){
  accessorFunction <- get(accessor)

  test_that(paste(accessor, "accepts multiple subjects and returns only those subjects"), {
    experiment <- metaclusteredSynthetic()
    requested <- c("subject1", "subject3")
    result <- suppressMessages(accessorFunction(experiment, subjects = requested))
    expect_setequal(unique(result$subject), requested)
  })

  test_that(paste(accessor, "accepts multiple cell subsets"), {
    experiment <- metaclusteredSynthetic()
    result <- suppressMessages(accessorFunction(experiment, cellSubsets = c("parent", "rareB")))
    expect_setequal(unique(result$cellSubset), c("parent", "rareB"))
  })

  test_that(paste(accessor, "accepts multiple markers and returns only those markers"), {
    experiment <- metaclusteredSynthetic()
    requested <- c("marker1", "marker2")
    result <- suppressMessages(accessorFunction(experiment, markers = requested))
    expect_true(all(requested %in% colnames(result)))
    otherMarkers <- setdiff(experiment$markerInfo$commonMarkerName, requested)
    expect_false(any(otherMarkers %in% colnames(result)))
  })

  test_that(paste(accessor, "accepts multiple metaclusters and returns only those"), {
    experiment <- metaclusteredSynthetic()
    result <- suppressMessages(accessorFunction(experiment, metaclusters = c(1, 2)))
    expect_true(all(result$metacluster %in% c(1, 2)))
    expect_gt(nrow(result), 0)
  })

  test_that(paste(accessor, "with single-value and default filters is unchanged"), {
    experiment <- metaclusteredSynthetic()
    allEvents <- suppressMessages(accessorFunction(experiment))
    expect_equal(nrow(allEvents), nrow(experiment$mergedExpr))
    oneSubject <- suppressMessages(accessorFunction(experiment, subjects = "subject2"))
    expect_identical(unique(oneSubject$subject), "subject2")
  })
}

test_that("makeMetaclusterHeatmaps accepts more than one child subset", {
  experiment <- metaclusteredSynthetic()
  outDir <- withr::local_tempdir()
  withr::local_dir(outDir)
  expect_no_error(suppressMessages(
    makeMetaclusterHeatmaps(experiment, childSubsets = c("rareA", "rareB"), filenamePrefix = "test", verbose = FALSE)
  ))
  expect_gt(length(list.files(outDir, pattern = "\\.png$")), 0)
})
