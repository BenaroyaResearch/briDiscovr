# Build a small discovrExperiment from synthetic FCS files written to a temporary
# directory. Each subject has a parent population (three well-separated
# phenotypes) and a rare subset that is a subset of the parent cells.
makeSyntheticExperiment <- function(nSubjects = 3, nCellsPerSubject = 600, nRareCells = 60,
                                    nMarkers = 6, seed = 1){
  set.seed(seed)
  fcsDir <- tempfile("syntheticFcs")
  dir.create(fcsDir)

  markerNames <- paste0("marker", seq_len(nMarkers))
  fcsMarkerNames <- paste0("Ch", seq_len(nMarkers), "_", markerNames)
  # three phenotypes: all markers low, medium, or high (on the arcsinh scale)
  phenotypeMeans <- c(1, 3, 5)

  writeFcs <- function(arcsinhValues, path){
    # undo the arcsinh(x * 0.2) transformation applied by setupDiscovrExperiment
    rawValues <- sinh(arcsinhValues) / 0.2
    colnames(rawValues) <- fcsMarkerNames
    ff <- flowCore::flowFrame(rawValues)
    params <- flowCore::parameters(ff)
    params$desc <- fcsMarkerNames
    flowCore::parameters(ff) <- params
    suppressWarnings(flowCore::write.FCS(ff, path))
  }

  fcsInfo <- data.frame()
  for(subjectIdx in seq_len(nSubjects)){
    subject <- paste0("subject", subjectIdx)
    phenotype <- sample(1:3, nCellsPerSubject, replace = TRUE)
    parentValues <- phenotypeMeans[phenotype] +
      matrix(rnorm(nCellsPerSubject * nMarkers, sd = 0.5), ncol = nMarkers)
    rareValues <- parentValues[sample(nCellsPerSubject, nRareCells), , drop = FALSE]

    parentPath <- file.path(fcsDir, paste0(subject, "_parent.fcs"))
    rarePath <- file.path(fcsDir, paste0(subject, "_rare.fcs"))
    writeFcs(parentValues, parentPath)
    writeFcs(rareValues, rarePath)
    fcsInfo <- rbind(
      fcsInfo,
      data.frame(subject = subject, cellSubset = c("parent", "rare"), filename = c(parentPath, rarePath))
    )
  }

  fcsInfoFile <- file.path(fcsDir, "fcsInfo.csv")
  write.csv(fcsInfo, fcsInfoFile, row.names = FALSE)
  markerInfoFile <- file.path(fcsDir, "markerInfo.csv")
  write.csv(
    data.frame(fixed = markerNames, desc = fcsMarkerNames, useToCluster = TRUE),
    markerInfoFile, row.names = FALSE
  )

  suppressMessages(setupDiscovrExperiment(
    markerInfoFile = markerInfoFile,
    fcsInfoFile = fcsInfoFile,
    parentPopulation = "parent",
    verbose = FALSE,
    checkMemory = FALSE
  ))
}
