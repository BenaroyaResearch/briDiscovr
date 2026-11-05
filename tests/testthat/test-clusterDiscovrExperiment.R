test_that("clusterDiscovrExperiment stores k and seed, and downstream steps keep them", {
  experiment <- makeSyntheticExperiment()
  expect_null(experiment$phenographClusterK)
  expect_null(experiment$clusterSeed)

  clustered <- suppressMessages(
    clusterDiscovrExperiment(experiment, k = 10, seed = 2024, verbose = FALSE)
  )
  expect_identical(clustered$phenographClusterK, 10)
  expect_identical(clustered$clusterSeed, 2024)
  expect_identical(clustered$clusterMethod, "phenograph")
  expect_output(print(clustered), "Clustering method: phenograph, k = 10, seed = 2024")

  normalized <- suppressMessages(normalizeDiscovrExperiment(clustered, verbose = FALSE))
  expect_identical(normalized$phenographClusterK, 10)
  expect_identical(normalized$clusterSeed, 2024)

  metaclustered <- suppressMessages(
    metaclusterDiscovrExperiment(normalized, nMetaclusters = 3, verbose = FALSE)
  )
  expect_identical(metaclustered$status, "metaclustered")
  expect_identical(metaclustered$phenographClusterK, 10)
  expect_identical(metaclustered$clusterSeed, 2024)

  recut <- suppressMessages(recutMetaclusters(metaclustered, nMetaclusters = 2))
  expect_identical(recut$phenographClusterK, 10)
  expect_identical(recut$clusterSeed, 2024)
})

test_that("the default seed is stored when seed is not given", {
  clustered <- suppressMessages(
    clusterDiscovrExperiment(makeSyntheticExperiment(), k = 10, verbose = FALSE)
  )
  expect_identical(clustered$clusterSeed, 12345)
})

test_that("re-clustering with different k and seed overwrites the stored values", {
  experiment <- makeSyntheticExperiment()
  first <- suppressMessages(clusterDiscovrExperiment(experiment, k = 10, seed = 1, verbose = FALSE))
  second <- suppressMessages(clusterDiscovrExperiment(experiment, k = 15, seed = 2, verbose = FALSE))
  expect_identical(c(first$phenographClusterK, first$clusterSeed), c(10, 1))
  expect_identical(c(second$phenographClusterK, second$clusterSeed), c(15, 2))
})

test_that("print handles clustered objects from older versions with no stored k or seed", {
  experiment <- suppressMessages(
    clusterDiscovrExperiment(makeSyntheticExperiment(), k = 10, verbose = FALSE)
  )
  experiment$phenographClusterK <- NULL
  experiment$clusterSeed <- NULL
  expect_output(print(experiment), "k = not recorded, seed = not recorded")
})

test_that("clusterDiscovrExperiment runs at small k with consecutive cluster IDs per sample", {
  # small k makes cells with no graph edges likely (previously an error or mislabeling)
  experiment <- makeSyntheticExperiment(nCellsPerSubject = 1500, seed = 2)
  clustered <- suppressMessages(clusterDiscovrExperiment(experiment, k = 5, verbose = FALSE))
  expect_false(any(is.na(clustered$mergedExpr$RPclust)))
  for(subject in unique(clustered$mergedExpr$samp)){
    labels <- clustered$mergedExpr$RPclust[clustered$mergedExpr$samp == subject]
    expect_identical(sort(unique(labels)), as.numeric(seq_len(max(labels))))
  }
})
