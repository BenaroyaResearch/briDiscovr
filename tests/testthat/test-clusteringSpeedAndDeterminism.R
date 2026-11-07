# Tests for the faster clustering steps (exact neighbor search with BiocNeighbors,
# rewritten jaccard_coeff, merged graph edges), the Leiden option, and seeding.

# plain-R version of the original Rphenograph jaccard_coeff(), used as a reference
referenceJaccard <- function(idx){
  nCols <- ncol(idx)
  weights <- matrix(0, nrow(idx) * nCols, 3)
  r <- 0
  for(i in seq_len(nrow(idx))){
    for(j in seq_len(nCols)){
      neighbor <- idx[i, j]
      shared <- length(intersect(idx[i, ], idx[neighbor, ]))
      if(shared > 0){
        r <- r + 1
        weights[r, ] <- c(i, neighbor, shared / (2 * nCols - shared) / 2)
      }
    }
  }
  weights
}

getNeighborsRann <- function(X, k){
  RANN::nn2(X, X, k + 1)[[1]][, -1]
}

# cluster IDs of each subject, as a named list
labelsBySubject <- function(experiment){
  split(experiment$mergedExpr$RPclust, experiment$mergedExpr$samp)
}

test_that("jaccard_coeff matches the original algorithm exactly", {
  set.seed(3); X <- matrix(rnorm(800 * 5), ncol = 5)
  for(k in c(3, 10)){
    neighbors <- getNeighborsRann(X, k)
    expect_identical(jaccard_coeff(neighbors), referenceJaccard(neighbors))
  }
})

test_that("BiocNeighbors finds the same neighbors as RANN, for any number of threads", {
  set.seed(4); X <- matrix(rnorm(3000 * 8), ncol = 8)
  oneThread <- briDiscovr:::findNearestNeighbors(X, k = 15, nThreads = 1)
  twoThreads <- briDiscovr:::findNearestNeighbors(X, k = 15, nThreads = 2)
  expect_identical(oneThread, twoThreads)
  expect_equal(oneThread, getNeighborsRann(X, 15), ignore_attr = TRUE)
})

test_that("merging duplicate edges keeps the same weighted connections", {
  set.seed(5); X <- matrix(rnorm(1500 * 6), ncol = 6)
  neighbors <- getNeighborsRann(X, 10)
  merged <- briDiscovr:::buildPhenographGraph(neighbors)

  # the unmerged graph, built as in briDiscovr 0.4.6 (but with vertices declared)
  links <- jaccard_coeff(neighbors)
  links <- links[links[, 1] > 0, ]
  unmerged <- igraph::graph_from_data_frame(
    data.frame(from = links[, 1], to = links[, 2], weight = links[, 3]),
    directed = FALSE, vertices = data.frame(name = seq_len(nrow(X)))
  )

  expect_lt(igraph::ecount(merged), igraph::ecount(unmerged))
  expect_equal(
    igraph::as_adjacency_matrix(merged, attr = "weight", sparse = FALSE),
    igraph::as_adjacency_matrix(unmerged, attr = "weight", sparse = FALSE)
  )
  partition <- sample(1:5, nrow(X), replace = TRUE)
  expect_equal(igraph::modularity(merged, partition), igraph::modularity(unmerged, partition))
})

test_that("clustering is identical across repeated runs and thread counts", {
  experiment <- makeSyntheticExperiment()
  for(method in c("phenograph", "leiden")){
    first <- suppressMessages(clusterDiscovrExperiment(
      experiment, method = method, k = 10, nThreads = 1, verbose = FALSE))
    second <- suppressMessages(clusterDiscovrExperiment(
      experiment, method = method, k = 10, nThreads = 2, verbose = FALSE))
    expect_identical(first$mergedExpr$RPclust, second$mergedExpr$RPclust)
  }
})

test_that("Leiden clustering runs and stores its settings; other methods clear them", {
  experiment <- makeSyntheticExperiment(nCellsPerSubject = 1500, seed = 2)
  leiden <- suppressMessages(clusterDiscovrExperiment(
    experiment, method = "leiden", k = 10, leidenResolution = 0.5, leidenIterations = 3, verbose = FALSE))
  expect_identical(leiden$clusterMethod, "leiden")
  expect_identical(leiden$leidenResolution, 0.5)
  expect_identical(leiden$leidenIterations, 3)
  expect_false(any(is.na(leiden$mergedExpr$RPclust)))
  for(labels in labelsBySubject(leiden)){
    expect_identical(sort(unique(labels)), as.numeric(seq_len(max(labels))))
  }
  expect_output(print(leiden), "Leiden resolution = 0.5, iterations = 3")

  # re-clustering with Louvain removes the Leiden settings
  louvain <- suppressMessages(clusterDiscovrExperiment(leiden, k = 10, verbose = FALSE))
  expect_null(louvain$leidenResolution)
  expect_null(louvain$leidenIterations)
})

test_that("unsupported methods and seeding options give errors", {
  experiment <- makeSyntheticExperiment()
  expect_error(clusterDiscovrExperiment(experiment, method = "kmeans", verbose = FALSE))
  expect_error(clusterDiscovrExperiment(experiment, setSeedOncePer = "subject", verbose = FALSE))
})

test_that("with setSeedOncePer = 'sample', a sample's clusters do not depend on other samples", {
  experiment <- makeSyntheticExperiment(nSubjects = 3, nCellsPerSubject = 1500, seed = 2)
  # the same experiment with the first subject removed
  withoutFirst <- experiment
  withoutFirst$mergedExpr <- experiment$mergedExpr[experiment$mergedExpr$samp != "subject1", ]

  for(method in c("phenograph", "leiden")){
    full <- suppressMessages(clusterDiscovrExperiment(
      experiment, method = method, k = 10, verbose = FALSE))
    reduced <- suppressMessages(clusterDiscovrExperiment(
      withoutFirst, method = method, k = 10, verbose = FALSE))
    expect_identical(full$clusterSetSeedOncePer, "sample")
    expect_identical(labelsBySubject(full)[c("subject2", "subject3")], labelsBySubject(reduced))
  }
})

test_that("setSeedOncePer = 'run' sets the seed once, as in earlier versions", {
  experiment <- makeSyntheticExperiment()
  byRun <- suppressMessages(clusterDiscovrExperiment(
    experiment, k = 10, setSeedOncePer = "run", verbose = FALSE))
  expect_identical(byRun$clusterSetSeedOncePer, "run")
  expect_output(print(byRun), "Seed set once per: run")

  # reproduce by hand: seed set once, then each subject clustered in turn
  set.seed(12345)
  expectedLabels <- lapply(split(experiment$mergedExpr, experiment$mergedExpr$samp)[unique(experiment$mergedExpr$samp)],
    function(subjectData){
      as.numeric(briDiscovr:::clusterEvents(
        subjectData[, experiment$clusteringMarkers], method = "phenograph", k = 10, nThreads = 1,
        leidenResolution = 1, leidenIterations = -1, verbose = FALSE))
    })
  expect_identical(labelsBySubject(byRun)[names(expectedLabels)], expectedLabels)
})
