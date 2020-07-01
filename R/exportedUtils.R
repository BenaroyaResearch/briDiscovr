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
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @importFrom utils View
#' @export
getSubjectCounts <- function(experiment){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object.",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'mergedExpr' %in% names(experiment)){
    stop(
      "The data passed to this function appears incorrectly formatted.",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'samp' %in% names(experiment$mergedExpr)){
    stop(
      "The data passed to this function appears incorrectly formatted.",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  # Check for donors with too many/too few collected events
  eventsByDonor <-
    experiment$mergedExpr %>%
    group_by(.data$samp) %>%
    summarise(n_events = n(), .groups = "drop_last") %>%
    as.data.frame()

  return(eventsByDonor)
}

#' Display the per-subject number of clusters in an experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment()}
#' @return A data frame containing the number of events detected for each subject
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise n
#' @importFrom utils View
#' @export
getSubjectClusters <- function(experiment){
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object.",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'mergedExpr' %in% names(experiment)){
    stop(
      "The data passed to this function appears incorrectly formatted.",
      "Please create your experiment using the 'setupDiscovrExperiment' function and try again."
    )
  }
  if (!'RPclust' %in% names(experiment$mergedExpr)){
    stop(
      "The data passed to this function does not appear to have been clustered yet.",
      "Please cluster your experiment using the 'clusterDiscovrExperiment' function and try again."
    )
  }
  # Get summary of the number of clusters generated for each subject
  nPhenoClusts <- experiment$mergedExpr %>%
    group_by(.data$samp) %>%
    summarize(k_clusters = max(.data$RPclust)) %>%
    as.data.frame()

  return(nPhenoClusts)
}
