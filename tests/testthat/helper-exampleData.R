# Helpers for tests that use the example dataset in tests/testthat/testdata/
# (10 donors, parent population "Tet-" and rare subset "Tet+"; see testdata/README.md).
#
# Running the full pipeline on this dataset takes a while (mainly the warpSet
# normalization), so each stage is computed once per test run and cached in
# `exampleDataCache`. Tests that modify an experiment get their own copy (R
# objects are copied on modification), so the cache is never changed.

exampleParentPopulation <- "Tet-"
exampleRareSubset <- "Tet+"
exampleSubjects <- paste0("donor", 1:10)
exampleClusterK <- 30
exampleClusterSeed <- 12345
exampleNMetaclusters <- 12

exampleDataCache <- new.env(parent = emptyenv())

# Evaluate an expression with its messages suppressed and its printed output
# discarded (warpSet normalization in flowStats prints progress with cat()).
quietly <- function(expr){
  utils::capture.output(result <- suppressMessages(expr))
  result
}

# Skip a test unless the example data and the packages needed for the full
# pipeline are available. The end-to-end tests take minutes, so skip on CRAN.
skipIfNoExampleData <- function(){
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("flowStats")
  testthat::skip_if_not(
    dir.exists(testthat::test_path("testdata", "fcs")),
    "example FCS files not found in tests/testthat/testdata/fcs"
  )
}

# Write a copy of fcsInfo.csv with absolute FCS paths, and return the paths of the
# input files to pass to setupDiscovrExperiment(). `subjects` keeps only those donors.
exampleDataPaths <- function(subjects = exampleSubjects, outDir = tempfile("exampleData")){
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  fcsInfo <- read.csv(testthat::test_path("testdata", "fcsInfo.csv"), stringsAsFactors = FALSE)
  fcsInfo <- fcsInfo[fcsInfo$subject %in% subjects, ]
  fcsInfo$filename <- normalizePath(testthat::test_path("testdata", "fcs", fcsInfo$filename), mustWork = TRUE)
  fcsInfoFile <- file.path(outDir, "fcsInfo.csv")
  write.csv(fcsInfo, fcsInfoFile, row.names = FALSE)
  list(
    fcsInfoFile = fcsInfoFile,
    markerInfoFile = normalizePath(testthat::test_path("testdata", "markerInfo.csv"), mustWork = TRUE),
    fcsInfo = fcsInfo
  )
}

# Run setupDiscovrExperiment() on the example data.
setupExampleExperiment <- function(subjects = exampleSubjects, ...){
  paths <- exampleDataPaths(subjects)
  suppressMessages(setupDiscovrExperiment(
    markerInfoFile = paths$markerInfoFile,
    fcsInfoFile = paths$fcsInfoFile,
    parentPopulation = exampleParentPopulation,
    verbose = FALSE,
    checkMemory = FALSE,
    ...
  ))
}

# Return the example experiment after the requested pipeline stage, computing
# (and caching) it and all earlier stages the first time it is requested.
# Stages: "initialized", "clustered", "normalized", "metaclustered".
getExampleExperiment <- function(stage = c("initialized", "clustered", "normalized", "metaclustered")){
  stage <- match.arg(stage)
  if(!is.null(exampleDataCache[[stage]])) return(exampleDataCache[[stage]])

  experiment <- switch(
    stage,
    initialized = setupExampleExperiment(),
    clustered = suppressMessages(clusterDiscovrExperiment(
      getExampleExperiment("initialized"),
      k = exampleClusterK, seed = exampleClusterSeed, nThreads = 1, verbose = FALSE
    )),
    normalized = quietly(normalizeDiscovrExperiment(
      getExampleExperiment("clustered"), verbose = FALSE
    )),
    metaclustered = suppressMessages(metaclusterDiscovrExperiment(
      getExampleExperiment("normalized"), nMetaclusters = exampleNMetaclusters, verbose = FALSE
    ))
  )
  exampleDataCache[[stage]] <- experiment
  experiment
}

# Read one example FCS file directly (untransformed), for comparison with package output.
readExampleFcs <- function(subject, cellSubset){
  fcsInfo <- exampleDataPaths(subject)$fcsInfo
  flowCore::read.FCS(
    fcsInfo$filename[fcsInfo$cellSubset == cellSubset],
    transformation = FALSE, truncate_max_range = FALSE
  )
}
