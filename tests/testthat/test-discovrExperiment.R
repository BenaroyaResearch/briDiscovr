# Tests for the discovrExperiment S3 class methods (R/discovrExperiment.R).
# Adapted from the pilot test suite in BenaroyaResearch/briDiscovr-2.
# Printing of clustering parameters is tested in test-clusterDiscovrExperiment.R.

test_that("is.discovrExperiment is TRUE only for discovrExperiment objects", {
  expect_true(is.discovrExperiment(create_mock_experiment()))
  expect_false(is.discovrExperiment(create_non_experiment()))
  expect_false(is.discovrExperiment(NULL))
  expect_false(is.discovrExperiment(data.frame(x = 1:5)))
  expect_false(is.discovrExperiment(1:10))
})

test_that("print.discovrExperiment shows markers, clustering markers, and status", {
  output <- paste(capture.output(print(create_mock_experiment(status = "initialized"))), collapse = "\n")
  expect_match(output, "An object of class 'discovrExperiment'", fixed = TRUE)
  expect_match(output, "All markers: CD4, CD8, CD3, CD19", fixed = TRUE)
  expect_match(output, "Clustering markers: CD4, CD8, CD3", fixed = TRUE)
  expect_match(output, "Experiment status: initialized", fixed = TRUE)
})

test_that("print.discovrExperiment omits clustering parameters before clustering", {
  output <- capture.output(print(create_mock_experiment()))
  expect_false(any(grepl("Clustering method", output)))
})

test_that("print.discovrExperiment errors for non-discovrExperiment input", {
  expect_error(
    print.discovrExperiment(create_non_experiment()),
    "must be an object of type 'discovrExperiment'"
  )
})

test_that("print.discovrExperiment warns about extra arguments", {
  expect_warning(
    capture.output(print(create_mock_experiment(), extra_arg = TRUE)),
    "Additional arguments were passed"
  )
})

test_that("getSampleNames returns each sample name once, in order of appearance", {
  expect_identical(getSampleNames(create_mock_experiment()), c("Sample1", "Sample2"))
})

test_that("getSampleNames errors for non-discovrExperiment input", {
  expect_error(
    getSampleNames(create_non_experiment()),
    "must be an object of type 'discovrExperiment'"
  )
})
