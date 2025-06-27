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
  # check inputs
  if(!is.data.frame(mergedExpr)) stop("Input 'mergedExpr' must be a data.frame")
  if (!all(sapply(mergedExpr[, !(colnames(mergedExpr) %in% groupCol)], is.numeric)))
    stop("Cannot perform warpSet normalization on non-numeric columns")
  
  # store information for later checking
  colnamesMergedExpr <- colnames(mergedExpr)
  groupColValuesMergedExpr <- mergedExpr[[groupCol]]
  
  # add a rowIdx to mergedExpr for easier merging later
  mergedExpr$rowIdx <- seq_len(nrow((mergedExpr)))
  
  # identify markers to normalize
  markersToNormalize <- setdiff(colnames(mergedExpr), c(groupCol, "rowIdx"))
  
  # normalize markers one at a time, and replace values in mergedExpr with these normalized values
  for(marker.tmp in markersToNormalize){
    normMergedExprByMarker <-
      normalizeWarpSetSingleMarker(
        mergedExpr[,c(marker.tmp, groupCol, "rowIdx")],
        marker = marker.tmp,
        peakNr = peakNr, groupCol = groupCol, seed = seed)
    mergedExpr[, marker.tmp] <-
      normMergedExprByMarker[
        match(mergedExpr$rowIdx, normMergedExprByMarker$rowIdx), marker.tmp]
  }
  
  # sort normalized data, check it, and remove rowIdx
  mergedExpr <- dplyr::arrange(mergedExpr, rowIdx)
  mergedExpr[["rowIdx"]] <- NULL
  if(!(identical(colnames(mergedExpr), colnamesMergedExpr) &
       identical(mergedExpr[[groupCol]], groupColValuesMergedExpr)))
    stop("Something went wrong with warpSet normalization")
  
  return(mergedExpr)
}

# Private function to apply warpSet normalization to a single marker
# this is done marker-by-marker to make the normalization more reproducible
# this function starts with a matrix containing sample information and one marker's expression
# it converts to flowSet splitting on a variable (generally "samp")
# should include only the marker to be normalized, the grouping column as groupCol, and rowIdx
# if peakNr (peak number) is specified for the current set of markers, it can be passed as peakNr
#' @importFrom flowStats warpSet
#' @importFrom rlang sym
normalizeWarpSetSingleMarker <-
  function(mergedExprByMarker, marker, peakNr = NULL, groupCol = "samp", seed = 12345){
    # check inputs
    if(!is.data.frame(mergedExprByMarker)) stop("Input 'mergedExpr' must be a data.frame")
    if(!setequal(colnames(mergedExprByMarker), c(marker, groupCol, "rowIdx")))
      stop("Input 'mergedExprByMarker' should contain only the marker, grouping column, and row index")
    
    # split mergedExpr by groupCol
    normMergedExprByMarker <- group_split(mergedExprByMarker, !!rlang::sym(groupCol), .keep = TRUE)
    
    # create flowSet
    flowSetMergedExprByMarker <-
      flowSet(
        lapply(normMergedExprByMarker,
               \(currFlowDat) flowFrame(as.matrix(currFlowDat[, !(colnames(currFlowDat) %in% c(groupCol, "rowIdx"))]))))
    
  # set random number seed before running warpSet, to ensure consistent results
  set.seed(seed)
  
  # run warpSet normalization
  flowSetMergedExprByMarker <-
    warpSet(flowSetMergedExprByMarker,
            stains = marker, peakNr = peakNr)
  
  # extract normalized expression from the flowSet object
  for(i in seq_len(length(normMergedExprByMarker))){
    normMergedExprByMarker[[i]][, marker] <-
      exprs(flowSetMergedExprByMarker[[i]])
  }
  
  # merge normalized data, check it
  normMergedExprByMarker <-
    bind_rows(normMergedExprByMarker) %>%
    dplyr::arrange(rowIdx)
  if(!(identical(colnames(mergedExprByMarker), colnames(normMergedExprByMarker)) &
       identical(mergedExprByMarker[[groupCol]], normMergedExprByMarker[[groupCol]])))
    stop("Something went wrong with warpSet normalization")
  
  return(normMergedExprByMarker)
}
