# Tests for metaclusterDiscovrExperiment() and recutMetaclusters()
# (R/metaclusteringPhase.R), using the example data.

# count events per subject, cell subset and metacluster directly from the events
countEventsByMetacluster <- function(experiment){
  events <- data.frame(
    subject = experiment$mergedExpr$samp,
    cellSubset = unlist(experiment$mergedExpr$cellSubset),
    metacluster = paste0("metacluster_", experiment$colIndices[experiment$mergedExpr$sampRpClust])
  )
  as.data.frame(table(events), stringsAsFactors = FALSE)
}

test_that("metaclustering sets status and assigns every retained cluster to one metacluster", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  expect_identical(metaclustered$status, "metaclustered")
  expect_identical(metaclustered$nMetaclusters, exampleNMetaclusters)
  expect_equal(metaclustered$kGroups, exampleNMetaclusters)
  expect_identical(metaclustered$linkage, "ward.D2")
  expect_identical(metaclustered$distance, "euclidean")
  expect_identical(metaclustered$metaclusterMarkers, metaclustered$clusteringMarkers)

  expect_false(anyDuplicated(names(metaclustered$colIndices)) > 0)
  expect_setequal(names(metaclustered$colIndices), unique(metaclustered$mergedExpr$sampRpClust))
  expect_setequal(unique(metaclustered$colIndices), seq_len(exampleNMetaclusters))
  # clustering settings from earlier steps are kept
  expect_identical(metaclustered$phenographClusterK, exampleClusterK)
})

test_that("clusters with fewer than pctInClusterThreshold percent of a subject's events are removed", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  metaclustered <- getExampleExperiment("metaclustered")
  clusterSizes <- table(paste0(clustered$mergedExpr$samp, "_", clustered$mergedExpr$RPclust))
  subjectSizes <- table(clustered$mergedExpr$samp)
  pctOfSubject <- 100 * as.vector(clusterSizes) / as.vector(subjectSizes[sub("_[0-9]+$", "", names(clusterSizes))])
  names(pctOfSubject) <- names(clusterSizes)

  expect_setequal(names(metaclustered$colIndices), names(pctOfSubject)[pctOfSubject >= 1])
  expect_equal(nrow(metaclustered$mergedExpr), sum(clusterSizes[pctOfSubject >= 1]))
  expect_equal(nrow(metaclustered$mergedExprNormalizedScaled), nrow(metaclustered$mergedExpr))
})

test_that("metacluster occupancy has every subject x metacluster combination", {
  skipIfNoExampleData()
  occupancy <- getExampleExperiment("metaclustered")$metaclusterOccupancy
  expect_named(occupancy, c("subject", "metacluster", "Tet-_frac", "Tet+_frac", "Tet-_cts", "Tet+_cts"))
  expect_equal(nrow(occupancy), length(exampleSubjects) * exampleNMetaclusters)
  expect_equal(nrow(unique(occupancy[, c("subject", "metacluster")])), nrow(occupancy))
})

test_that("metacluster occupancy counts match the events in each metacluster", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  occupancy <- metaclustered$metaclusterOccupancy
  counts <- countEventsByMetacluster(metaclustered)
  for(subset in c(exampleParentPopulation, exampleRareSubset)){
    expected <- counts[counts$cellSubset == subset, ]
    observed <- merge(occupancy, expected, by = c("subject", "metacluster"))
    expect_equal(nrow(observed), nrow(occupancy))
    expect_equal(observed[[paste0(subset, "_cts")]], observed$Freq)
  }
})

test_that("metacluster occupancy fractions sum to 1 within each subject and cell subset", {
  skipIfNoExampleData()
  occupancy <- getExampleExperiment("metaclustered")$metaclusterOccupancy
  for(subset in c(exampleParentPopulation, exampleRareSubset)){
    fracColumn <- paste0(subset, "_frac")
    ctsColumn <- paste0(subset, "_cts")
    sums <- tapply(occupancy[[fracColumn]], occupancy$subject, sum)
    expect_equal(unname(as.vector(sums)), rep(1, length(sums)), tolerance = 1e-8)
    expect_true(all(occupancy[[fracColumn]] >= 0 & occupancy[[fracColumn]] <= 1))
    # each fraction is the count divided by the subject's total for that subset
    totals <- tapply(occupancy[[ctsColumn]], occupancy$subject, sum)
    expect_equal(occupancy[[fracColumn]], occupancy[[ctsColumn]] / as.vector(totals[occupancy$subject]))
  }
})

test_that("getMetaclusterOccupancy returns the stored occupancy, rounded", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  result <- getMetaclusterOccupancy(metaclustered, precision = 3)
  expect_equal(result$`Tet+_frac`, round(metaclustered$metaclusterOccupancy$`Tet+_frac`, 3))
  expect_identical(result$`Tet+_cts`, metaclustered$metaclusterOccupancy$`Tet+_cts`)
})

test_that("dropSamples removes those samples everywhere", {
  skipIfNoExampleData()
  dropped <- suppressMessages(metaclusterDiscovrExperiment(
    getExampleExperiment("normalized"), dropSamples = c("donor3", "donor7"), nMetaclusters = 8, verbose = FALSE
  ))
  for(element in c("mergedExpr", "mergedExprNormalizedScaled", "clusterMeans", "clusterMeansNormalizedScaled")){
    expect_false(any(c("donor3", "donor7") %in% dropped[[element]]$samp), info = element)
  }
  expect_false(any(c("donor3", "donor7") %in% dropped$metaclusterOccupancy$subject))
  expect_false(any(grepl("^donor[37]_", names(dropped$colIndices))))
  expect_equal(nrow(dropped$metaclusterOccupancy), 8 * 8)
})

test_that("dropMarkers removes markers from metaclustering only", {
  skipIfNoExampleData()
  dropped <- suppressMessages(metaclusterDiscovrExperiment(
    getExampleExperiment("normalized"), dropMarkers = c("CD38", "IgM"), verbose = FALSE
  ))
  expect_identical(dropped$metaclusterMarkers, setdiff(dropped$clusteringMarkers, c("CD38", "IgM")))
  expect_true(all(c("CD38", "IgM") %in% rownames(dropped$allSubsetAllSubjectNormalizedScaled)))
})

test_that("a stricter pctInClusterThreshold keeps fewer clusters", {
  skipIfNoExampleData()
  strict <- suppressMessages(metaclusterDiscovrExperiment(
    getExampleExperiment("normalized"), pctInClusterThreshold = 5, verbose = FALSE
  ))
  expect_identical(strict$pctInClusterThreshold, 5)
  expect_lt(length(strict$colIndices), length(getExampleExperiment("metaclustered")$colIndices))
  expect_true(all(names(strict$colIndices) %in% names(getExampleExperiment("metaclustered")$colIndices)))
})

test_that("metaclustering a clustered experiment normalizes it first, with a message", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  # use zScore everywhere so the automatic normalization is fast
  clustered$markerInfo$normalizationMethod <- "zScore"
  expect_message(
    result <- metaclusterDiscovrExperiment(clustered, nMetaclusters = 5, verbose = FALSE),
    "Normalization will now be performed prior to metaclustering"
  )
  expect_identical(result$status, "metaclustered")
  expect_false(is.null(result$mergedExprNormalizedScaled))
})

test_that("metaclustering supports other linkage methods", {
  skipIfNoExampleData()
  averageLinkage <- suppressMessages(metaclusterDiscovrExperiment(
    getExampleExperiment("normalized"), linkage = "average", nMetaclusters = 6, verbose = FALSE
  ))
  expect_identical(averageLinkage$linkage, "average")
  expect_equal(averageLinkage$kGroups, 6)
  # average linkage builds a different tree from the default (Ward's method)
  expect_false(identical(averageLinkage$clusterDendrogram$merge, getExampleExperiment("metaclustered")$clusterDendrogram$merge))
})

test_that("metaclusterDiscovrExperiment validates its input", {
  expect_error(metaclusterDiscovrExperiment(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(
    metaclusterDiscovrExperiment(create_mock_experiment(status = "initialized")),
    "must have status 'clustered' or 'normalized'"
  )
})

# recutMetaclusters -----------------------------------------------------------

test_that("recutting at the same number of metaclusters reproduces the original results", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  recut <- suppressMessages(recutMetaclusters(metaclustered, nMetaclusters = exampleNMetaclusters))
  expect_identical(recut$colIndices, metaclustered$colIndices)
  expect_equal(recut$metaclusterOccupancy, metaclustered$metaclusterOccupancy)
})

test_that("recutting cuts the same dendrogram at a different height", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  recut <- suppressMessages(recutMetaclusters(metaclustered, nMetaclusters = 5))
  expect_identical(recut$nMetaclusters, 5)
  expect_equal(recut$kGroups, 5)
  expect_equal(recut$clusterDendrogram$merge, metaclustered$clusterDendrogram$merge)
  expect_identical(names(recut$colIndices), names(metaclustered$colIndices))
  # each of the 12 original metaclusters lies entirely within one of the 5 new ones
  nesting <- unique(data.frame(original = metaclustered$colIndices, recut = recut$colIndices))
  expect_equal(nrow(nesting), exampleNMetaclusters)
  # occupancy is recomputed for the new metaclusters, and counts are unchanged in total
  expect_equal(nrow(recut$metaclusterOccupancy), length(exampleSubjects) * 5)
  expect_equal(sum(recut$metaclusterOccupancy$`Tet+_cts`), sum(metaclustered$metaclusterOccupancy$`Tet+_cts`))
  counts <- countEventsByMetacluster(recut)
  observed <- merge(recut$metaclusterOccupancy, counts[counts$cellSubset == exampleRareSubset, ])
  expect_equal(observed$`Tet+_cts`, observed$Freq)
})

test_that("recutMetaclusters validates its input", {
  expect_error(recutMetaclusters(create_non_experiment()), "not a valid DISCOV-R experiment object")
  expect_error(
    recutMetaclusters(create_mock_experiment(status = "normalized")),
    "must have already been metaclustered"
  )
  missingDendrogram <- create_mock_experiment(status = "metaclustered")
  expect_error(recutMetaclusters(missingDendrogram), "does not have the data expected")
})
