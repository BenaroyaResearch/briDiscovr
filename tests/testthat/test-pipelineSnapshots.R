# Snapshot (regression) tests of the full pipeline on the example data.
#
# These record key results the first time they run (in _snaps/) and fail if later
# runs give different results. They check that results do not change silently;
# they do not check that the results are correct. Other test files check that.
#
# Clustering results depend on the igraph version, and normalized values (and so
# metaclusters) on the flowStats version, so snapshots are stored separately for
# each combination of their major.minor versions. When a change in results is
# expected (e.g. an intentional change to clustering), review the differences
# and accept them with testthat::snapshot_review() or testthat::snapshot_accept().

snapshotVariant <- function(){
  majorMinor <- function(package) paste(unlist(utils::packageVersion(package))[1:2], collapse = ".")
  paste0("igraph-", majorMinor("igraph"), "_flowStats-", majorMinor("flowStats"))
}

test_that("number of clusters per subject is unchanged", {
  skipIfNoExampleData()
  expect_snapshot_value(
    getSubjectClusters(getExampleExperiment("clustered")),
    style = "json2", variant = snapshotVariant()
  )
})

test_that("cluster sizes are unchanged", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  clusterSizes <- as.data.frame(table(sample = clustered$mergedExpr$samp, cluster = clustered$mergedExpr$RPclust))
  clusterSizes <- clusterSizes[clusterSizes$Freq > 0, ]
  clusterSizes <- clusterSizes[order(clusterSizes$sample, as.numeric(as.character(clusterSizes$cluster))), ]
  rownames(clusterSizes) <- NULL
  expect_snapshot_value(
    lapply(clusterSizes, as.character),
    style = "json2", variant = snapshotVariant()
  )
})

test_that("normalized expression distributions are unchanged", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  quantiles <- sapply(
    normalized$mergedExprNormalizedScaled[, sort(normalized$markersNormalized)],
    stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
  )
  expect_snapshot_value(round(quantiles, 3), style = "json2", variant = snapshotVariant(), tolerance = 1e-3)
})

test_that("metacluster assignments of clusters are unchanged", {
  skipIfNoExampleData()
  colIndices <- getExampleExperiment("metaclustered")$colIndices
  expect_snapshot_value(
    as.list(colIndices[order(names(colIndices))]),
    style = "json2", variant = snapshotVariant()
  )
})

test_that("metacluster occupancy is unchanged", {
  skipIfNoExampleData()
  occupancy <- getMetaclusterOccupancy(getExampleExperiment("metaclustered"), precision = 4)
  occupancy <- occupancy[order(occupancy$subject, occupancy$metacluster), ]
  rownames(occupancy) <- NULL
  expect_snapshot_value(occupancy, style = "json2", variant = snapshotVariant(), tolerance = 1e-4)
})
