# Regression tests for PhenoGraph graph construction when some cells have no
# Jaccard edges of their own (small k). Previously these cells could cause an
# error (cell missing from the graph) or silently misalign cluster labels.

getNeighbors <- function(X, k){
  RANN::nn2(X, X, k + 1)[[1]][, -1]
}

test_that("graph has one vertex per cell, in row order, when some cells appear nowhere", {
  # at k = 5 this dataset has cells with no edge at all (previously the error case)
  set.seed(1); X <- matrix(rnorm(2000 * 9), ncol = 9)
  g <- briDiscovr:::buildPhenographGraph(getNeighbors(X, k = 5))
  expect_equal(igraph::vcount(g), nrow(X))
  expect_identical(as.integer(igraph::V(g)$name), seq_len(nrow(X)))
  expect_gt(sum(igraph::degree(g, loops = FALSE) == 0), 0)
})

test_that("graph has one vertex per cell, in row order, when edgeless cells appear only as targets", {
  # at k = 7 this dataset previously gave a full-size but misordered graph (silent case)
  set.seed(2); X <- matrix(rnorm(1000 * 9), ncol = 9)
  g <- briDiscovr:::buildPhenographGraph(getNeighbors(X, k = 7))
  expect_equal(igraph::vcount(g), nrow(X))
  expect_identical(as.integer(igraph::V(g)$name), seq_len(nrow(X)))
})

test_that("isolated cells get the majority label of connected neighbors, and IDs are consecutive", {
  set.seed(1); X <- matrix(rnorm(2000 * 9), ncol = 9)
  neighbors <- getNeighbors(X, k = 5)
  g <- briDiscovr:::buildPhenographGraph(neighbors)
  set.seed(12345)
  louvainLabels <- as.integer(igraph::membership(igraph::cluster_louvain(g)))
  result <- briDiscovr:::assignIsolatedCells(louvainLabels, g, neighbors)

  isolated <- which(igraph::degree(g, loops = FALSE) == 0)
  expect_equal(result$nIsolated, length(isolated))
  expect_length(result$membership, nrow(X))
  expect_identical(sort(unique(result$membership)), seq_len(max(result$membership)))

  # connected cells keep the same partition (labels may be renumbered):
  # each old label maps to exactly one new label and vice versa
  connected <- setdiff(seq_len(nrow(X)), isolated)
  labelPairs <- unique(data.frame(old = louvainLabels[connected], new = result$membership[connected]))
  expect_equal(nrow(labelPairs), length(unique(labelPairs$old)))
  expect_equal(nrow(labelPairs), length(unique(labelPairs$new)))

  # each isolated cell with connected neighbors has a label that is the most common among them
  for(cell in isolated){
    connectedNeighbors <- neighbors[cell, ][!(neighbors[cell, ] %in% isolated)]
    if(length(connectedNeighbors) == 0) next
    neighborCounts <- table(result$membership[connectedNeighbors])
    expect_true(neighborCounts[as.character(result$membership[cell])] == max(neighborCounts))
  }
})

test_that("assignIsolatedCells breaks ties toward the closest neighbor and handles all-isolated neighbors", {
  # 5 cells; cell 1 is isolated; neighbors listed closest first
  g <- igraph::graph_from_data_frame(
    data.frame(from = c(2, 4), to = c(3, 5)),
    directed = FALSE,
    vertices = data.frame(name = 1:5)
  )
  neighbors <- rbind(
    c(4, 2), # cell 1: closest is cell 4 (cluster 2), then cell 2 (cluster 1): tie -> cluster 2
    c(3, 1),
    c(2, 1),
    c(5, 1),
    c(4, 1)
  )
  labels <- c(3L, 1L, 1L, 2L, 2L)
  result <- briDiscovr:::assignIsolatedCells(labels, g, neighbors)
  expect_identical(result$membership, c(2L, 1L, 1L, 2L, 2L))
  expect_equal(result$nIsolated, 1)

  # an isolated cell whose neighbors are all isolated keeps its own cluster
  gAllIsolated <- igraph::make_empty_graph(3, directed = FALSE)
  resultAllIsolated <- briDiscovr:::assignIsolatedCells(
    c(1L, 2L, 3L), gAllIsolated, rbind(c(2, 3), c(1, 3), c(1, 2))
  )
  expect_identical(resultAllIsolated$membership, c(1L, 2L, 3L))
})

test_that("Rphenograph membership lines up with input rows at small k", {
  set.seed(2); X <- matrix(rnorm(1000 * 9), ncol = 9)
  set.seed(12345)
  out <- suppressMessages(capture.output(res <- Rphenograph(X, k = 7)))
  expect_length(igraph::membership(res[[2]]), nrow(X))
  expect_identical(as.integer(igraph::V(res[[1]])$name), seq_len(nrow(X)))

  set.seed(1); X5 <- matrix(rnorm(2000 * 9), ncol = 9)
  set.seed(12345)
  out <- suppressMessages(capture.output(res5 <- Rphenograph(X5, k = 5)))
  expect_length(igraph::membership(res5[[2]]), nrow(X5))
})
