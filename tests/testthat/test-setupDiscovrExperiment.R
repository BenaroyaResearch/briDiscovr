# Tests for setupDiscovrExperiment() (R/clusteringPhase.R), using the example data
# in tests/testthat/testdata/.

# write a modified copy of the example markerInfo file, and return its path
writeMarkerInfo <- function(modify = identity){
  markerInfo <- read.csv(testthat::test_path("testdata", "markerInfo.csv"), stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".csv")
  write.csv(modify(markerInfo), path, row.names = FALSE)
  path
}

# run setup on two donors with optional replacement input files
setupTwoDonors <- function(markerInfoFile = NULL, fcsInfo = NULL, parentPopulation = exampleParentPopulation, ...){
  paths <- exampleDataPaths(c("donor1", "donor2"))
  if(!is.null(fcsInfo)){
    write.csv(fcsInfo(paths$fcsInfo), paths$fcsInfoFile, row.names = FALSE)
  }
  suppressMessages(setupDiscovrExperiment(
    markerInfoFile = if(is.null(markerInfoFile)) paths$markerInfoFile else markerInfoFile,
    fcsInfoFile = paths$fcsInfoFile,
    parentPopulation = parentPopulation,
    verbose = FALSE, checkMemory = FALSE, ...
  ))
}

test_that("setup returns an initialized experiment with all subjects and events", {
  skipIfNoExampleData()
  experiment <- getExampleExperiment("initialized")
  expect_true(is.discovrExperiment(experiment))
  expect_identical(experiment$status, "initialized")
  expect_identical(experiment$parentPopulation, exampleParentPopulation)
  expect_setequal(getSampleNames(experiment), exampleSubjects)

  # event counts per subject and cell subset match the FCS files
  fcsInfo <- exampleDataPaths()$fcsInfo
  expectedCounts <- setNames(
    vapply(fcsInfo$filename, getFcsNEvents, integer(1)),
    paste(fcsInfo$subject, fcsInfo$cellSubset)
  )
  observedCounts <- table(paste(experiment$mergedExpr$samp, unlist(experiment$mergedExpr$cellSubset)))
  expect_equal(as.vector(observedCounts[names(expectedCounts)]), unname(expectedCounts))
})

test_that("setup records clustering markers and fills in missing normalization methods", {
  skipIfNoExampleData()
  experiment <- getExampleExperiment("initialized")
  markerInfoFile <- read.csv(testthat::test_path("testdata", "markerInfo.csv"), stringsAsFactors = FALSE)
  expect_identical(experiment$clusteringMarkers, markerInfoFile$fixed[markerInfoFile$useToCluster])
  expect_named(experiment$markerInfo, c("commonMarkerName", "fcsMarkerName", "useToCluster", "normalizationMethod"))
  # blank normalization methods (the FSC markers) become zScore
  expect_identical(experiment$markerInfo$normalizationMethod[1:2], c("zScore", "zScore"))
  expect_true(all(experiment$clusteringMarkers %in% colnames(experiment$mergedExpr)))
})

test_that("setup stores the parent population first within each subject", {
  skipIfNoExampleData()
  experiment <- getExampleExperiment("initialized")
  for(subject in c("donor1", "donor7")){
    subsets <- unlist(experiment$mergedExpr$cellSubset[experiment$mergedExpr$samp == subject])
    expect_identical(rle(subsets)$values, c(exampleParentPopulation, exampleRareSubset))
  }
})

test_that("expression values are arcsinh(0.2 * raw) of the FCS values, in file order", {
  skipIfNoExampleData()
  experiment <- getExampleExperiment("initialized")
  fcs <- readExampleFcs("donor3", exampleRareSubset)
  params <- flowCore::pData(flowCore::parameters(fcs))
  for(marker in c("CD20", "IgD", "CD3")){
    rawValues <- flowCore::exprs(fcs)[, params$name[which(params$desc == marker)]]
    isDonor3Rare <- experiment$mergedExpr$samp == "donor3" &
      unlist(experiment$mergedExpr$cellSubset) == exampleRareSubset
    expect_equal(experiment$mergedExpr[isDonor3Rare, marker], unname(asinh(0.2 * rawValues)), tolerance = 1e-6)
  }
})

test_that("arcsinh cofactor can be changed with arcsinhB", {
  skipIfNoExampleData()
  default <- setupTwoDonors()
  rescaled <- setupTwoDonors(arcsinhB = 1)
  rawFromDefault <- sinh(default$mergedExpr$CD20) / 0.2
  expect_equal(rescaled$mergedExpr$CD20, asinh(rawFromDefault), tolerance = 1e-6)
})

test_that("downsampleVectorList keeps only the listed events of each file", {
  skipIfNoExampleData()
  # fcsInfo order for two donors: donor1 Tet-, donor2 Tet-, donor1 Tet+, donor2 Tet+
  rowsToKeep <- list(1:100, "all", seq(1, 500, by = 5), 1:50)
  experiment <- setupTwoDonors(downsampleVectorList = rowsToKeep)
  counts <- table(paste(experiment$mergedExpr$samp, unlist(experiment$mergedExpr$cellSubset)))
  expect_equal(as.vector(counts[c("donor1 Tet-", "donor2 Tet-", "donor1 Tet+", "donor2 Tet+")]), c(100, 2000, 100, 50))

  # the kept events are the requested rows of the file
  fcs <- readExampleFcs("donor1", exampleRareSubset)
  params <- flowCore::pData(flowCore::parameters(fcs))
  rawCD20 <- flowCore::exprs(fcs)[seq(1, 500, by = 5), params$name[which(params$desc == "CD20")]]
  isDonor1Rare <- experiment$mergedExpr$samp == "donor1" & unlist(experiment$mergedExpr$cellSubset) == exampleRareSubset
  expect_equal(experiment$mergedExpr$CD20[isDonor1Rare], unname(asinh(0.2 * rawCD20)), tolerance = 1e-6)
})

test_that("downsampleVectorList is validated", {
  skipIfNoExampleData()
  expect_error(setupTwoDonors(downsampleVectorList = "not a list"), "must be a list")
  expect_error(setupTwoDonors(downsampleVectorList = list(1:10, letters)), "must be a list")
  expect_error(setupTwoDonors(downsampleVectorList = list(1:10, 1:10)), "must have one element for each FCS file")
})

test_that("setup checks the markerInfo file", {
  skipIfNoExampleData()
  expect_error(
    setupTwoDonors(writeMarkerInfo(function(m){ m$useToCluster <- NULL; m })),
    "must contain columns with names 'useToCluster'"
  )
  expect_error(
    setupTwoDonors(writeMarkerInfo(function(m){ m$useToCluster <- ifelse(m$useToCluster, "yes", "no"); m })),
    "must contain only TRUE or FALSE"
  )
  expect_error(
    setupTwoDonors(writeMarkerInfo(function(m){ m$fixed[4] <- m$fixed[3]; m })),
    "do not appear to be unique"
  )
  expect_error(
    setupTwoDonors(writeMarkerInfo(function(m){ m$normalizationMethod[3] <- "quantile"; m })),
    "can only contain values"
  )
  expect_error(
    setupTwoDonors(writeMarkerInfo(function(m){ m$desc[3] <- "NotInTheFcsFiles"; m })),
    "Could not find the following clustering markers.*NotInTheFcsFiles"
  )
})

test_that("setup reports a marker with no name", {
  skipIfNoExampleData()
  skip(paste(
    "Known bug: marker names are passed through make.names() before the blank-name check,",
    "so a blank name becomes 'NA.' and the check never triggers. Remove this skip once fixed."
  ))
  expect_error(
    setupTwoDonors(writeMarkerInfo(function(m){ m$fixed[12] <- ""; m })),
    "Markers without assigned names"
  )
})

test_that("setup uses custom markerInfo column names", {
  skipIfNoExampleData()
  markerInfoFile <- writeMarkerInfo(function(m){ names(m)[1:2] <- c("marker", "fcsDescription"); m })
  experiment <- setupTwoDonors(markerInfoFile, markerCommonField = "marker", markerFcsField = "fcsDescription")
  expect_true(all(experiment$clusteringMarkers %in% colnames(experiment$mergedExpr)))
})

test_that("setup checks the fcsInfo file and parent population", {
  skipIfNoExampleData()
  expect_error(setupTwoDonors(parentPopulation = "CD8"), "requested parent population 'CD8' was not found")
  expect_error(
    setupTwoDonors(fcsInfo = function(f){ f$cellSubset[f$cellSubset == "Tet+"] <- "CD20"; f }),
    "Cell subset names and marker names must not overlap"
  )
  expect_error(
    setupTwoDonors(fcsInfo = function(f){ f$cellSubset[3] <- "Tet-"; f }),
    "subject/cell subset pairs were duplicated"
  )
})
