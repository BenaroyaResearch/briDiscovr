## Copyright (C) 2025  Matt Dufort, Mario Rosasco, Virginia Muir, and Benaroya Research Institute
##
## This file is part of the briDiscovr package

# Private function to read in data from FCS files and associate the data with
# subject IDs. Accepts a data frame with at least two columns:'filename',
# which should contain paths to .fcs files, and the field indicated by
# 'indexField' which will be used as the index in the returned list
buildFcsList <- function(fileFrame, indexField = "subject", ...) {
  stopifnot(all(c(indexField, "filename") %in% colnames(fileFrame)))
  stopifnot(is.data.frame(fileFrame))
  tmpList <- list()
  for (currIdx in unique(fileFrame[[indexField]])) {
    currFcsName = as.character(fileFrame[which(fileFrame[[indexField]]==currIdx), "filename"])
    if(file.info(currFcsName)$size > 0){
      tmpList[[currIdx]] <- flowCore::read.FCS(currFcsName, transformation = FALSE, ...)
    }
  }
  return(tmpList)
}

# Private function to process FCSList objects to ensure that all data include all required markers with correct names
processFcsList <- function(fcs, markerInfo){
  # confirm that all markers to be used in clustering are mappable to fcs parameter names
  clusteringMarkerDesc <- markerInfo[markerInfo$useToCluster, "fcsMarkerName", drop = TRUE]
  missingMarkers <- clusteringMarkerDesc[!clusteringMarkerDesc %in% pData(parameters(fcs))$desc]
  
  # check for missingMarkers in name field of flowFrame parameters, where desc is empty
  missingMarkersToRename <- intersect(missingMarkers, pData(parameters(fcs))$name)
  missingMarkersToRename <-
    missingMarkersToRename[
      is.na(pData(parameters(fcs))$desc[match(missingMarkersToRename, pData(parameters(fcs))$name)])]
  if (length(missingMarkersToRename) > 0){
    pData(parameters(fcs))$desc[match(missingMarkers, pData(parameters(fcs))$name)] <-
      missingMarkersToRename
  }
  missingMarkers <- setdiff(missingMarkers, missingMarkersToRename)
  
  if(length(missingMarkers) > 0){
    stop("Could not find the following clustering markers in the .fcs data\n", paste0(missingMarkers, collapse = ", "))
  }
  
  # Tidy marker names (does this cause markers without a value in pData(parameters(fcs))$desc to be dropped?)
  pData(parameters(fcs))$desc <-
    markerInfo$commonMarkerName[match(pData(parameters(fcs))$desc, markerInfo$fcsMarkerName)]
  
  # This changes parameters(fcs)$name, featureNames(fcs), and colnames(fcs) - aka events colnames - all in one fell swoop.
  # note colnames has to be the one from flowCore
  flowCore::colnames(fcs) <- make.names(pData(parameters(fcs))$desc)
  
  # Remove markers that aren't informative/shared between panels (i.e. duplicated NAs)
  fcs <- fcs[,!(duplicated(flowCore::colnames(fcs)) | duplicated(flowCore::colnames(fcs), fromLast = TRUE))]
  fcs <- fcs[, order(flowCore::colnames(fcs))]
  
  return(fcs)
}

# Private function to generate a list of colors of arbitrary length.
getColorList <- function(n, seed = 42){
  set.seed(seed)
  qualColorPals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  allColors = unlist(mapply(RColorBrewer::brewer.pal, qualColorPals$maxcolors, rownames(qualColorPals)))

  retList = c()

  if (n > length(allColors)){
    message("Note: ", n, " colors requested, ", length(allColors), " unique colors available. Some colors will be repeated.")
    while(n > length(allColors)){
      retList = c(retList, sample(allColors))
      n = n - length(allColors)
    }
  }
  retList = c(retList, sample(allColors)[1:n])
  return(retList)
}

# Private function to check normalizationMethod within markerInfo against a list of possibilities
# if normalizationMethod is not found, set it to "zScore"
checkMarkerInfoNormalizationMethod <- function(markerInfo){
  # check normalizationMethod, and set to zScore if not present
  if("normalizationMethod" %in% names(markerInfo)){
    if(
      !all(
        (markerInfo$normalizationMethod %in% c("zScore", "none", "warpSet")) |
        str_detect(markerInfo$normalizationMethod, "^warpSet[0-9]+$") |
        is.na(markerInfo$normalizationMethod))){
      stop("Column 'normalizationMethod' in markerInfo file can only contain values 'zScore', 'none', 'warpset[#]', or blanks/NAs. Please check the contents of your markerInfo file.")
    }
    markerInfo$normalizationMethod[is.na(markerInfo$normalizationMethod)] <- "zScore"
  } else {
    markerInfo$normalizationMethod <- "zScore"
  }
  return(markerInfo)
}

# Private function to apply warpSet normalization to a set of markers from mergedExpr
# this function starts with a matrix, converts to flowSet splitting on a variable (generally "samp")
# should include only the markers to be normalized and the grouping column as groupCol
# if peakNr (peak number) is specified for the current set of markers, it can be passed as peakNr
#' @importFrom flowStats warpSet
#' @importFrom rlang sym
normalizeWarpSetMergedExpr <- function(mergedExpr, peakNr = NULL, groupCol = "samp", seed = 12345){
  if (!all(sapply(mergedExpr[, !(colnames(mergedExpr) %in% groupCol)], is.numeric)))
    stop("Cannot perform warpSet normalization on non-numeric columns")
  
  # add a rowIdx to mergedExpr for easier merging later
  mergedExpr$rowIdx <- seq_len(nrow((mergedExpr)))
  
  # split mergedExpr by groupCol
  normMergedExpr <- group_split(mergedExpr, !!rlang::sym(groupCol), .keep = TRUE)
  
  # create flowSet
  flowSetMergedExpr <-
    flowSet(lapply(normMergedExpr, \(currFlowDat) flowFrame(as.matrix(currFlowDat[, !(colnames(currFlowDat) %in% c(groupCol, "rowIdx"))]))))
  
  # set random number seed before running warpSet, to ensure consistent results
  set.seed(seed)
  
  # run warpSet normalization
  flowSetMergedExpr <-
    warpSet(flowSetMergedExpr, stains = setdiff(colnames(mergedExpr), c(groupCol, "rowIdx")), peakNr = peakNr)
  
  # extract normalized expression from the flowSet object
  for(i in seq_len(length(normMergedExpr))){
    normMergedExpr[[i]][, match(colnames(exprs(flowSetMergedExpr[[i]])), colnames(normMergedExpr[[i]]))] <-
      exprs(flowSetMergedExpr[[i]])
  }
  
  # merge normalized data, check it, and remove rowIdx
  normMergedExpr <-
    bind_rows(normMergedExpr) %>%
    dplyr::arrange(rowIdx)
  if(!(identical(colnames(mergedExpr), colnames(normMergedExpr)) &
       identical(mergedExpr[[groupCol]], normMergedExpr[[groupCol]])))
    stop("Something went wrong with warpSet normalization")
  normMergedExpr[["rowIdx"]] <- NULL
  
  return(normMergedExpr)
}
