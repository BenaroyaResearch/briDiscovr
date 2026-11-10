# Tests for runUmapDiscovrExperiment() and plotDensityNormalizedExprsDiscovrExperiment()
# (R/exportedUtils.R), using the example data. Input validation with mock objects
# is in test-exportedUtils.R.

# count retained events per sample and cell subset in a data frame with samp and cellSubset columns
countBySampleAndSubset <- function(data){
  table(paste(data$samp, unlist(data$cellSubset)))
}

test_that("UMAP by frequency keeps 1 of every n events per sample and cell subset", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  umapData <- suppressMessages(runUmapDiscovrExperiment(
    metaclustered, downsampleFreq = c(parentPopulation = 20, childPopulations = 10), seed = 1
  ))
  expect_s3_class(umapData, "data.frame")
  expect_named(umapData, c("samp", "cellSubset", "metacluster", "UMAP1", "UMAP2"))
  available <- countBySampleAndSubset(metaclustered$mergedExprNormalizedScaled)
  freq <- ifelse(grepl(exampleRareSubset, names(available), fixed = TRUE), 10, 20)
  expect_equal(as.vector(countBySampleAndSubset(umapData)[names(available)]), as.vector(ceiling(available / freq)))
  expect_true(all(umapData$metacluster %in% seq_len(exampleNMetaclusters)))
  expect_false(anyNA(umapData[, c("UMAP1", "UMAP2")]))
})

test_that("UMAP by number keeps at most n events per sample and cell subset", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  umapData <- suppressMessages(runUmapDiscovrExperiment(
    metaclustered, downsampleBy = "number",
    downsampleNumber = c(parentPopulation = 60, childPopulations = 30), seed = 1
  ))
  counts <- countBySampleAndSubset(umapData)
  expect_true(all(counts[grepl(exampleParentPopulation, names(counts), fixed = TRUE)] == 60))
  expect_true(all(counts[grepl(exampleRareSubset, names(counts), fixed = TRUE)] == 30))
})

test_that("UMAP is reproducible with a seed", {
  skipIfNoExampleData()
  metaclustered <- getExampleExperiment("metaclustered")
  runSmallUmap <- function(seed){
    suppressMessages(runUmapDiscovrExperiment(
      metaclustered, downsampleFreq = c(parentPopulation = 50, childPopulations = 40), seed = seed
    ))
  }
  first <- runSmallUmap(seed = 5)
  expect_identical(runSmallUmap(seed = 5), first)
  expect_false(identical(runSmallUmap(seed = 6)$UMAP1, first$UMAP1))
})

test_that("UMAP options return the UMAP object, the expression values, and chosen markers", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  result <- suppressMessages(runUmapDiscovrExperiment(
    normalized, umapMarkers = c("CD20", "CD27", "IgD"),
    downsampleFreq = c(parentPopulation = 50, childPopulations = 40),
    seed = 1, returnUmapObject = TRUE, returnExpressionNormalizedScaled = TRUE
  ))
  expect_named(result, c("data", "umapObject"))
  expect_s3_class(result$umapObject, "umap")
  expect_true(all(c("CD20", "CD27", "IgD", "UMAP1", "UMAP2") %in% colnames(result$data)))
  # a normalized (not metaclustered) experiment has no metacluster column
  expect_false("metacluster" %in% colnames(result$data))
  expect_equal(nrow(result$umapObject$layout), nrow(result$data))
})

test_that("UMAP normalizes a clustered experiment first, using markerInfo", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  clustered$markerInfo$normalizationMethod <- "zScore"
  expect_message(
    result <- runUmapDiscovrExperiment(
      clustered, downsampleFreq = c(parentPopulation = 50, childPopulations = 50), seed = 1
    ),
    "has not yet had marker expression normalized"
  )
  # nothing was removed by metaclustering, so all events of the clustered experiment are available
  expected <- ceiling(countBySampleAndSubset(clustered$mergedExpr) / 50)
  expect_equal(nrow(result), sum(expected))
})

# plotDensityNormalizedExprsDiscovrExperiment --------------------------------

test_that("density plots are returned for each clustering marker, or all markers", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  plots <- plotDensityNormalizedExprsDiscovrExperiment(normalized)
  expect_named(plots, normalized$clusteringMarkers)
  expect_true(all(vapply(plots, ggplot2::is.ggplot, logical(1))))
  expect_match(plots$IgD$labels$title, "IgD (normalization: warpSet3)", fixed = TRUE)

  allPlots <- plotDensityNormalizedExprsDiscovrExperiment(normalized, clusteringMarkersOnly = FALSE)
  expect_setequal(names(allPlots), normalized$markersNormalized)
})

test_that("density plots can be written to a PDF file", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  outFile <- file.path(withr::local_tempdir(), "densities.pdf")
  plotDensityNormalizedExprsDiscovrExperiment(normalized, filenameOut = outFile)
  expect_true(file.exists(outFile))
  expect_identical(readBin(outFile, "raw", 4), charToRaw("%PDF"))
})

test_that("density plots can be written to PNG, JPEG, and TIFF files", {
  skipIfNoExampleData()
  normalized <- getExampleExperiment("normalized")
  outDir <- withr::local_tempdir()
  # first bytes of each file format
  signatures <- list(
    png = as.raw(c(0x89, 0x50, 0x4e, 0x47)),
    jpeg = as.raw(c(0xff, 0xd8, 0xff)),
    tiff = NULL # "II*" or "MM*", depending on byte order
  )
  for(format in names(signatures)){
    outFile <- file.path(outDir, paste0("densities.", format))
    plotDensityNormalizedExprsDiscovrExperiment(normalized, filenameOut = outFile, width = 4, height = 2, res = 50)
    expect_true(file.exists(outFile), info = format)
    if(!is.null(signatures[[format]])){
      expect_identical(readBin(outFile, "raw", length(signatures[[format]])), signatures[[format]], info = format)
    }
  }
})

test_that("density plots normalize a clustered experiment first", {
  skipIfNoExampleData()
  clustered <- getExampleExperiment("clustered")
  clustered$markerInfo$normalizationMethod <- "zScore"
  expect_message(
    plots <- plotDensityNormalizedExprsDiscovrExperiment(clustered, verbose = FALSE),
    "Normalization will be run prior to plotting"
  )
  expect_named(plots, clustered$clusteringMarkers)
})
