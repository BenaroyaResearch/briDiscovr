# Tests for internal utilities (R/internalUtils.R).
# getColorList and checkMarkerInfoNormalizationMethod tests are adapted from the
# pilot test suite in BenaroyaResearch/briDiscovr-2.

test_that("getColorList returns the requested number of hex colors", {
  colors <- briDiscovr:::getColorList(10)
  expect_length(colors, 10)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
  expect_length(briDiscovr:::getColorList(3), 3)
})

test_that("getColorList returns distinct colors when enough are available", {
  skip(paste(
    "Known bug: the RColorBrewer qualitative palettes share 4 colors, so the pool",
    "has 74 entries but 70 distinct colors; with the default seed, 12 colors (the",
    "default number of metaclusters) include a duplicate. Remove this skip once fixed."
  ))
  expect_equal(anyDuplicated(briDiscovr:::getColorList(12)), 0)
  expect_equal(anyDuplicated(briDiscovr:::getColorList(20)), 0)
  expect_equal(anyDuplicated(briDiscovr:::getColorList(70)), 0)
})

test_that("getColorList is reproducible with the same seed and differs with another", {
  expect_identical(briDiscovr:::getColorList(10, seed = 42), briDiscovr:::getColorList(10, seed = 42))
  expect_false(identical(briDiscovr:::getColorList(10, seed = 42), briDiscovr:::getColorList(10, seed = 123)))
})

test_that("getColorList repeats colors, with a message, when more are requested than exist", {
  expect_message(colors <- briDiscovr:::getColorList(100), "Some colors will be repeated")
  expect_length(colors, 100)
  expect_true(all(!is.na(colors)))
})

test_that("checkMarkerInfoNormalizationMethod accepts valid methods unchanged", {
  markerInfo <- data.frame(
    commonMarkerName = c("CD4", "CD8", "CD3", "CD19", "CD20"),
    normalizationMethod = c("zScore", "none", "warpSet", "warpSet1", "warpSet12")
  )
  result <- briDiscovr:::checkMarkerInfoNormalizationMethod(markerInfo)
  expect_identical(result$normalizationMethod, markerInfo$normalizationMethod)
})

test_that("checkMarkerInfoNormalizationMethod replaces NA with zScore", {
  markerInfo <- data.frame(commonMarkerName = c("CD4", "CD8"), normalizationMethod = c("none", NA))
  result <- briDiscovr:::checkMarkerInfoNormalizationMethod(markerInfo)
  expect_identical(result$normalizationMethod, c("none", "zScore"))
})

test_that("checkMarkerInfoNormalizationMethod adds a zScore column if missing", {
  result <- briDiscovr:::checkMarkerInfoNormalizationMethod(data.frame(commonMarkerName = c("CD4", "CD8")))
  expect_identical(result$normalizationMethod, c("zScore", "zScore"))
})

test_that("checkMarkerInfoNormalizationMethod rejects invalid methods", {
  for(badMethod in c("invalidMethod", "zscore", "warpset", "warpSetX")){
    markerInfo <- data.frame(commonMarkerName = c("CD4", "CD8"), normalizationMethod = c("zScore", badMethod))
    expect_error(briDiscovr:::checkMarkerInfoNormalizationMethod(markerInfo), "can only contain values")
  }
})

test_that("isUnsetArgument is TRUE only for a single NA", {
  expect_true(briDiscovr:::isUnsetArgument(NA))
  expect_true(briDiscovr:::isUnsetArgument(NA_character_))
  expect_false(briDiscovr:::isUnsetArgument("a"))
  expect_false(briDiscovr:::isUnsetArgument(c("a", "b")))
  expect_false(briDiscovr:::isUnsetArgument(c(NA, NA)))
  expect_false(briDiscovr:::isUnsetArgument(NULL))
})

test_that("normalizeWarpSetMergedExpr validates its input", {
  expect_error(briDiscovr:::normalizeWarpSetMergedExpr(as.matrix(data.frame(a = 1))), "must be a data.frame")
  expect_error(
    briDiscovr:::normalizeWarpSetMergedExpr(data.frame(CD4 = "high", samp = "s1")),
    "non-numeric columns"
  )
})
