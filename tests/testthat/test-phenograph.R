# Tests for Rphenograph() and find_neighbors() (R/phenograph.R).
# Input-validation and find_neighbors tests are adapted from the pilot test suite
# in BenaroyaResearch/briDiscovr-2. The pilot skipped full-clustering tests because
# of label misalignment with small datasets; that was fixed in briDiscovr 0.4.7
# (see test-phenographGraph.R), so those tests are included here.

test_that("Rphenograph validates its input", {
  data <- matrix(rnorm(100), ncol = 4)
  expect_error(Rphenograph("not a matrix", k = 5), "Wrong input data")
  expect_error(Rphenograph(data, k = 0), "k must be a positive integer")
  expect_error(Rphenograph(data, k = -5), "k must be a positive integer")
  expect_error(Rphenograph(matrix(rnorm(40), ncol = 4), k = 15), "k must be smaller than the total number of points")
})

test_that("Rphenograph returns a graph and communities with one entry per row", {
  X <- withr::with_seed(42, rbind(
    matrix(rnorm(150 * 4, mean = 0), ncol = 4),
    matrix(rnorm(150 * 4, mean = 6), ncol = 4)
  ))
  result <- withr::with_seed(12345, suppressMessages(capture.output(res <- Rphenograph(X, k = 15))))
  expect_length(res, 2)
  expect_s3_class(res[[1]], "igraph")
  expect_equal(igraph::vcount(res[[1]]), nrow(X))
  expect_length(igraph::membership(res[[2]]), nrow(X))
})

test_that("Rphenograph separates two well-separated groups", {
  X <- withr::with_seed(42, rbind(
    matrix(rnorm(150 * 4, mean = 0), ncol = 4),
    matrix(rnorm(150 * 4, mean = 6), ncol = 4)
  ))
  withr::with_seed(12345, suppressMessages(capture.output(res <- Rphenograph(X, k = 15))))
  membership <- as.integer(igraph::membership(res[[2]]))
  # no community contains points from both groups
  communitiesGroup1 <- unique(membership[1:150])
  communitiesGroup2 <- unique(membership[151:300])
  expect_length(intersect(communitiesGroup1, communitiesGroup2), 0)
})

test_that("Rphenograph accepts a data.frame and gives the same result as a matrix", {
  X <- withr::with_seed(42, matrix(rnorm(200 * 4), ncol = 4))
  withr::with_seed(12345, suppressMessages(capture.output(fromMatrix <- Rphenograph(X, k = 10))))
  withr::with_seed(12345, suppressMessages(capture.output(fromDataFrame <- Rphenograph(as.data.frame(X), k = 10))))
  expect_identical(
    as.integer(igraph::membership(fromMatrix[[2]])),
    as.integer(igraph::membership(fromDataFrame[[2]]))
  )
})

test_that("Rphenograph is reproducible with the same seed", {
  X <- withr::with_seed(42, matrix(rnorm(300 * 5), ncol = 5))
  withr::with_seed(12345, suppressMessages(capture.output(first <- Rphenograph(X, k = 10))))
  withr::with_seed(12345, suppressMessages(capture.output(second <- Rphenograph(X, k = 10))))
  expect_identical(
    as.integer(igraph::membership(first[[2]])),
    as.integer(igraph::membership(second[[2]]))
  )
})

test_that("find_neighbors returns a k-column matrix of valid row indices, self first", {
  X <- withr::with_seed(42, matrix(rnorm(100), ncol = 4))
  neighbors <- find_neighbors(X, k = 5)
  expect_true(is.matrix(neighbors))
  expect_equal(dim(neighbors), c(25, 5))
  expect_equal(neighbors[, 1], 1:25)
  expect_true(all(neighbors >= 1 & neighbors <= 25))
})

test_that("find_neighbors orders neighbors from closest to farthest", {
  X <- withr::with_seed(42, matrix(rnorm(200), ncol = 4))
  neighbors <- find_neighbors(X, k = 6)
  distances <- as.matrix(dist(X))
  for(i in seq_len(nrow(X))){
    neighborDistances <- distances[i, neighbors[i, ]]
    expect_false(is.unsorted(neighborDistances))
    # no point outside the neighbor set is closer than the farthest neighbor
    expect_true(all(distances[i, -neighbors[i, ]] >= max(neighborDistances)))
  }
})
