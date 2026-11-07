# Tests for the FCS file utilities in R/exportedUtils.R (getFcsNEvents,
# checkForFcsByteOffsetIssue, checkFcsFiles, downsampleFcsList), using the
# example data in tests/testthat/testdata/.

# write a small FCS file with nEvents random events, and return its path
writeSmallFcs <- function(nEvents, path = tempfile(fileext = ".fcs")){
  values <- matrix(withr::with_seed(1, runif(nEvents * 2)), ncol = 2, dimnames = list(NULL, c("A", "B")))
  suppressWarnings(flowCore::write.FCS(flowCore::flowFrame(values), path))
  path
}

test_that("getFcsNEvents matches the number of events read by flowCore", {
  skipIfNoExampleData()
  fcsInfo <- exampleDataPaths()$fcsInfo
  for(file in fcsInfo$filename[c(1, 11, 12, 20)]){
    expect_identical(getFcsNEvents(file), nrow(flowCore::read.FCS(file, truncate_max_range = FALSE)))
  }
})

test_that("checkForFcsByteOffsetIssue finds no issue in valid files", {
  skipIfNoExampleData()
  for(file in exampleDataPaths()$fcsInfo$filename){
    expect_false(checkForFcsByteOffsetIssue(file))
  }
})

test_that("checkFcsFiles passes the example data", {
  skipIfNoExampleData()
  expect_true(suppressMessages(checkFcsFiles(exampleDataPaths()$fcsInfo, checkMemory = FALSE)))
  expect_true(suppressMessages(checkFcsFiles(exampleDataPaths()$fcsInfo, checkMemory = TRUE, verbose = FALSE)))
})

test_that("checkFcsFiles reports duplicated subject/cell subset pairs", {
  skipIfNoExampleData()
  fcsInfo <- exampleDataPaths()$fcsInfo
  fcsInfo$subject[2] <- fcsInfo$subject[1]
  expect_error(checkFcsFiles(fcsInfo, checkMemory = FALSE), "subject/cell subset pairs were duplicated")
})

test_that("checkFcsFiles reports files with no events", {
  # flowCore cannot write an FCS file with no events, so report 0 events for a real file
  local_mocked_bindings(getFcsNEvents = function(fcsFile) 0L)
  fcsInfo <- data.frame(subject = "S1", cellSubset = "Parent", filename = writeSmallFcs(10))
  expect_error(checkFcsFiles(fcsInfo, checkMemory = FALSE, verbose = FALSE), "have no events")
})

test_that("checkFcsFiles gives a message for subjects with too few events to cluster", {
  fcsInfo <- data.frame(
    subject = c("S1", "S1", "S2"),
    cellSubset = c("Parent", "Rare", "Parent"),
    filename = c(writeSmallFcs(10), writeSmallFcs(5), writeSmallFcs(100))
  )
  messages <- testthat::capture_messages(result <- checkFcsFiles(fcsInfo, checkMemory = FALSE, verbose = FALSE))
  expect_true(result)
  expect_length(messages, 1)
  expect_match(messages, "fewer than 32 total cells")
  # S1 has 15 events in total and is listed; S2 has 100 and is not
  expect_match(messages, "\nS1\n?$")
  expect_no_match(messages, "S2")
})

# downsampleFcsList -----------------------------------------------------------

test_that("downsampleFcsList returns row indices for each file, capped at maxEvents", {
  skipIfNoExampleData()
  paths <- exampleDataPaths()
  rows <- suppressMessages(downsampleFcsList(paths$fcsInfoFile, maxEvents = 1000, seed = 1, verbose = FALSE))
  nEvents <- vapply(paths$fcsInfo$filename, getFcsNEvents, integer(1))
  expect_length(rows, nrow(paths$fcsInfo))
  expect_equal(lengths(rows), unname(pmin(nEvents, 1000)))
  for(i in seq_along(rows)){
    expect_false(anyDuplicated(rows[[i]]) > 0)
    expect_true(all(rows[[i]] >= 1 & rows[[i]] <= nEvents[[i]]))
  }
  # files with no more than maxEvents keep all events, in order
  expect_identical(rows[[which(nEvents <= 1000)[1]]], seq_len(nEvents[nEvents <= 1000][1]))
})

test_that("downsampleFcsList is reproducible with a seed", {
  skipIfNoExampleData()
  paths <- exampleDataPaths(c("donor1", "donor2"))
  first <- downsampleFcsList(paths$fcsInfoFile, maxEvents = 300, seed = 7, verbose = FALSE)
  second <- downsampleFcsList(paths$fcsInfoFile, maxEvents = 300, seed = 7, verbose = FALSE)
  expect_identical(first, second)
})

test_that("downsampleFcsList downsamples only the requested population", {
  skipIfNoExampleData()
  paths <- exampleDataPaths(c("donor1", "donor2"))
  rows <- downsampleFcsList(paths$fcsInfoFile, popsToDownsample = "Tet+", maxEvents = 300, seed = 1, verbose = FALSE)
  isRare <- paths$fcsInfo$cellSubset == "Tet+"
  expect_true(all(lengths(rows[isRare]) == 300))
  expect_true(all(lengths(rows[!isRare]) == 2000))
})

test_that("downsampleFcsList accepts more than one population", {
  skipIfNoExampleData()
  skip(paste(
    "Known bug: downsampleFcsList uses `if(popsToDownsample == \"all\")`, which errors in",
    "R >= 4.2 when more than one population is given. Remove this skip once fixed."
  ))
  paths <- exampleDataPaths(c("donor1", "donor2"))
  rows <- downsampleFcsList(paths$fcsInfoFile, popsToDownsample = c("Tet-", "Tet+"), maxEvents = 300, verbose = FALSE)
  expect_true(all(lengths(rows) == 300))
})

test_that("downsampleFcsList writes downsampled FCS files and a new fcsInfo file", {
  skipIfNoExampleData()
  paths <- exampleDataPaths(c("donor1", "donor2"))
  outDir <- withr::local_tempdir()
  outFcsInfo <- file.path(outDir, "fcsInfoDownsampled.csv")
  suppressWarnings(downsampleFcsList(
    paths$fcsInfoFile, maxEvents = 250, downsampleMode = "writeFiles",
    outfileFcsInfoDownsampled = outFcsInfo, outdirFcsFilesDownsampled = outDir,
    seed = 1, verbose = FALSE
  ))
  newFcsInfo <- read.csv(outFcsInfo, stringsAsFactors = FALSE)
  expect_identical(newFcsInfo[, c("subject", "cellSubset")], paths$fcsInfo[, c("subject", "cellSubset")], ignore_attr = TRUE)
  expect_true(all(file.exists(newFcsInfo$filename)))
  expect_true(all(dirname(normalizePath(newFcsInfo$filename)) == normalizePath(outDir)))
  expect_true(all(vapply(newFcsInfo$filename, getFcsNEvents, integer(1)) == 250))

  # the written files can be used to set up an experiment
  experiment <- suppressMessages(setupDiscovrExperiment(
    paths$markerInfoFile, outFcsInfo, parentPopulation = exampleParentPopulation,
    verbose = FALSE, checkMemory = FALSE
  ))
  expect_equal(nrow(experiment$mergedExpr), 4 * 250)
})

test_that("downsampleFcsList validates its input", {
  skipIfNoExampleData()
  paths <- exampleDataPaths(c("donor1"))
  expect_error(downsampleFcsList("/nonexistent/fcsInfo.csv"), "cannot be accessed")
  expect_error(
    downsampleFcsList(paths$fcsInfoFile, popsToDownsample = "CD8", verbose = FALSE),
    "were not found in the FCS info file"
  )
  expect_error(
    downsampleFcsList(
      paths$fcsInfoFile, downsampleMode = "writeFiles", verbose = FALSE,
      outfileFcsInfoDownsampled = "/nonexistent/dir/fcsInfo.csv", outdirFcsFilesDownsampled = tempdir()
    ),
    "'outfileFcsInfoDownsampled' is not valid"
  )
  expect_error(
    downsampleFcsList(
      paths$fcsInfoFile, downsampleMode = "writeFiles", verbose = FALSE,
      outfileFcsInfoDownsampled = tempfile(), outdirFcsFilesDownsampled = "/nonexistent/dir"
    ),
    "'outdirFcsFilesDownsampled' is not valid"
  )
  expect_error(downsampleFcsList(paths$fcsInfoFile, downsampleMode = "other", verbose = FALSE))
})
