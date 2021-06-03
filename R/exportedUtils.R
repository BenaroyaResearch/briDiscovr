## Copyright (C) 2020  Mario Rosasco, Virginia Muir, and Benaroya Research Institute
##
## This file is part of the briDiscovr package


#' Check for byte offset corruption in the header of an .fcs file
#'
#' Opens an .fcs file in binary mode and reads the 2 values of
#' the byte offset from the text header. These should agree with one another as
#' a check for file data integrity, but in practice the values are often off
#' by one. This seems to be related to exporting rare gated populations from
#' flowJo, and can often be corrected by re-gating and re-exporting.
#'
#' @param fcsFile A string containing the location of an .fcs file
#' @return Boolean: TRUE if there's an offset issue detected, FALSE otherwise
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom stringr str_extract str_replace str_replace_all regex
#' @export
checkForFcsByteOffsetIssue <- function(fcsFile){
  fileHandle <- file(fcsFile, "rb")
  # toss un-needed header info
  readChar(fileHandle, nchars = 26)

  headDataStartOffset <- as.integer(readChar(fileHandle, 8))
  headDataEndOffset <- as.integer(readChar(fileHandle, 8))

  # check for start/end validity
  if(!(headDataEndOffset > headDataStartOffset) | (headDataStartOffset < 42)){
    close(fileHandle)
    return(TRUE)
  }

  tmpData <- readChar(fileHandle, headDataStartOffset - 42)
  close(fileHandle)

  # get the separator character
  sepChar <- str_extract(tmpData, regex("\\$BEGINDATA.", dotall = T)) %>% str_replace("\\$BEGINDATA", "")

  textDataStartOffset <-
    str_extract(tmpData, paste0("\\$BEGINDATA\\", sepChar, "[0-9]+")) %>%
    str_replace_all(paste0("\\$BEGINDATA|\\", sepChar), "") %>%
    as.integer()
  textDataEndOffset <-
    str_extract(tmpData, paste0("\\$ENDDATA\\", sepChar, "[0-9]+")) %>%
    str_replace_all(paste0("\\$ENDDATA\\", sepChar), "") %>%
    as.integer()

  byteOffsetOk <-
    (textDataStartOffset == headDataStartOffset) &
    (textDataEndOffset == headDataEndOffset)

  return(!byteOffsetOk)
}

#' Check the number of events in .fcs files
#'
#' Many clustering algorithms cannot handle data objects with no events, so this
#' utility is provided as a way to quickly check the number of events in a file.
#' The number of events is parsed from the text header using low-level file
#' I/O, so counting the events in very large .fcs files will not take
#' longer than smaller files.
#'
#' @param fcsFile A string containing the location of an .fcs file
#' @return Numeric indicating the number of events listed in the .fcs header
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom stringr str_extract str_replace str_replace_all regex
#' @export
getFcsNEvents <- function(fcsFile) {
  fileHandle <- file(fcsFile, "rb")
  # toss un-needed header info
  readChar(fileHandle, nchars = 26)

  headDataStartOffset <- as.integer(readChar(fileHandle, 8))
  headDataEndOffset <- as.integer(readChar(fileHandle, 8))
  tmpData <- readChar(fileHandle, headDataStartOffset - 42)
  close(fileHandle)

  # get the separator character
  sepChar <- str_extract(tmpData, regex("\\$TOT.", dotall = T)) %>% str_replace("\\$TOT", "")

  nEvents <-
    str_extract(tmpData, paste0("\\$TOT\\", sepChar, "[0-9]+")) %>%
    str_replace_all(paste0("\\$TOT|\\", sepChar), "") %>%
    as.integer()

  return(nEvents)
}


#' Display the per-subject event counts in an experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment()}
#' @return A data frame containing the number of events detected for each subject
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n rename
#' @export
getSubjectCounts <- function(experiment){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'mergedExpr' %in% names(experiment)){
    stop(
      "The data passed to this function appears incorrectly formatted. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'samp' %in% names(experiment$mergedExpr)){
    stop(
      "The data passed to this function appears incorrectly formatted. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  # Check for donors with too many/too few collected events
  eventsByDonor <-
    experiment$mergedExpr %>%
    rename(sample = .data$samp) %>%
    group_by(.data$sample) %>%
    summarise(nEvents = n(), .groups = "drop_last") %>%
    as.data.frame()

  return(eventsByDonor)
}

#' Display the per-subject number of clusters in an experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment()}
#' @return A data frame containing the number of events detected for each subject
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @export
getSubjectClusters <- function(experiment){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'mergedExpr' %in% names(experiment)){
    stop(
      "The data passed to this function appears incorrectly formatted. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'RPclust' %in% names(experiment$mergedExpr)){
    stop(
      "The data passed to this function does not appear to have been clustered yet. ",
      "Please cluster your experiment using the 'clusterDiscovrExperiment' function and try again."
    )
  }
  # Get summary of the number of clusters generated for each subject
  nPhenoClusts <- experiment$mergedExpr %>%
    rename(sample = .data$samp) %>%
    group_by(.data$sample) %>%
    summarize(kClusters = max(.data$RPclust)) %>%
    as.data.frame()

  return(nPhenoClusts)
}

#' Display the fraction of events that are found in each metacluster for each
#' sample and cell subset
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment()}
#' @param precision An integer indicating the number of significant digits to return for fractional values
#' @return A data frame containing the fraction and number of events in each metacluster
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom stringr str_detect
#' @export
getMetaclusterOccupancy <- function(experiment, precision=2){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'metaclusterOccupancy' %in% names(experiment)){
    stop(
      "The data passed to this function does not appear to have been metaclustered yet. ",
      "Please perform metaclustering using the 'metaclusterDiscovrExperiment' function and try again."
    )
  }
  if (!is.numeric(precision) || precision %% 1 != 0){
    stop("The 'precision' argument must be set to an integer number.")
  }

  metaxVals = experiment$metaclusterOccupancy
  fracCols = stringr::str_detect(colnames(metaxVals), "_frac")
  metaxVals[,fracCols] <- sapply(metaxVals[,fracCols], round, digits = precision)

  return(as.data.frame(metaxVals))
}

#' Compute and return individual event z-scores from a metaclustered experiment
#'
#' @param experiment A metaclustered discovrExperiment
#' @param cellSubsets A string or vector of strings indicating which cell subsets to return event z-scores from
#' @param metaclusters An integer or vector of integers indicating which metaclusters to return event z-scores from
#' @param subjects A string or vector of strings indicating which subjects to return event z-scores from
#' @param markers A string or vector of strings indicating which markers to return event z-scores from
#' @return A data frame containing event z-scores (computed per-subject) from the events satisfying the requested filters. For any filter field without a specified value, all possible events will be returned.
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
getEventZScores <- function(experiment, cellSubsets = NA, metaclusters = NA, subjects = NA, markers = NA){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'mergedExpr' %in% names(experiment)){
    stop(
      "The data passed to this function appears incorrectly formatted. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'RPclust' %in% names(experiment$mergedExpr)){
    stop(
      "The data passed to this function does not appear to have been clustered yet. ",
      "Please cluster your experiment using the 'clusterDiscovrExperiment' function and try again."
    )
  }
  if(experiment$status != "metaclustered"){
    stop(
      "The current experiment status is ", experiment$status,
      ". The experiment must be metaclustered before using this function."
    )
  }

  # If arguments are left as NA assign to all values
  if(is.na(cellSubsets)){
    cellSubsets <- unique(experiment$mergedExpr$cellSubset)
  }
  if(is.na(subjects)){
    subjects <- unique(experiment$mergedExpr$samp)
  }
  if(is.na(markers)){
    markers <- unique(experiment$markerInfo$commonMarkerName)
  }
  if(is.na(metaclusters)){
    metaclusters <- 1:experiment$nMetaclusters
  }

  # Check validity of filters
  if(!all(cellSubsets %in% experiment$mergedExpr$cellSubset)){
    stop(
      "Requested values in the 'cellSubset' filter could not be found in the data. ",
      "Please check these values and try again."
    )
  }
  if(!all(subjects %in% experiment$mergedExpr$samp)){
    stop(
      "Requested values in the 'subject' filter could not be found in the data. ",
      "Please check these values and try again."
    )
  }
  if(!all(markers %in% names(experiment$mergedExpr))){
    stop(
      "Requested values in the 'markers' filter could not be found in the data. ",
      "Please check these values and try again."
    )
  }
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5){ abs(x - round(x)) < tol }
  if(!all(is.wholenumber(metaclusters))){
    stop(
      "Requested values in the 'metaclusters' filter do not appear to be whole number values. ",
      "Please check these values and try again."
    )
  }
  if(max(metaclusters) > experiment$nMetaclusters){
    stop(
      "Requested values in the 'metaclusters' filter are greater than the number of metaclusters (", experiment$nMetaclusters, "). ",
      "Please check these values and try again."
    )
  }

  message("Computing Z-scores for the selected events. Please be patient, this may take a moment...")

  # get std dev and mean in long format for merge with event data
  message("Computing per-subject standard deviations...")
  parentStdDev <-
    experiment$mergedExpr %>%
    dplyr::select(-.data$cellSubset, -.data$RPclust, -.data$sampRpClust) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(.data$samp) %>%
    dplyr::summarise_all(sd) %>%
    tidyr::gather("marker", "subjectStdDev", -samp)

  message("Computing per-subject means...")
  parentMean <-
    experiment$mergedExpr %>%
    dplyr::select(-.data$cellSubset, -.data$RPclust, -.data$sampRpClust) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(.data$samp) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather("marker", "subjectMean", -samp)

  parentMeanStdDev <-
    inner_join(parentMean, parentStdDev, by = c("samp", "marker"))

  # get metacluster indices
  metaclusterIndices <-
    data.frame(metacluster = experiment$colIndices) %>%
    tibble::rownames_to_column("sampRpClust") %>%
    dplyr::filter(.data$metacluster %in% metaclusters)

  # filter the data and make into long format to compute z-scores
  message("Filtering data and computing event z-scores...")
  rowsToUse <- (experiment$mergedExpr$samp %in% subjects) & (experiment$mergedExpr$cellSubset %in% cellSubsets)
  colsToUse <- unique(c(markers, "samp", "cellSubset", "sampRpClust"))

  zScoreData <-
    experiment$mergedExpr[rowsToUse, colsToUse] %>%
    # track event ID to help re-spread table later
    dplyr::mutate(event = row_number()) %>%
    tidyr::gather("marker", "value", -.data$samp, -.data$cellSubset, -.data$sampRpClust, -.data$event) %>%
    # merge in metacluster indices
    dplyr::inner_join(metaclusterIndices, by = "sampRpClust") %>%
    # merge in stdev and mean
    dplyr::inner_join(parentMeanStdDev, by = c("samp", "marker")) %>%
    # compute z-score
    dplyr::mutate(zScore = (.data$value-.data$subjectMean)/.data$subjectStdDev) %>%
    # reformat data
    dplyr::select(.data$samp, .data$marker, .data$zScore, .data$metacluster, .data$event) %>%
    tidyr::spread("marker", "zScore") %>%
    dplyr::rename(subject = .data$samp) %>%
    dplyr::select(-.data$event)

  return(zScoreData)
}
