## Copyright (C) 2022  Matt Dufort, Mario Rosasco, Virginia Muir, and Benaroya Research Institute
##
## This file is part of the briDiscovr package


#' Check for byte offset corruption in the header of an .fcs file
#'
#' Uses flowCore functions to check the 2 values of the byte offset from the
#' text header. These should agree with one another as a check for file data
#' integrity, but in practice the values are often off by one. This seems to
#' be related to exporting rare gated populations from flowJo, and can often
#' be corrected by re-gating and re-exporting. This function previously used
#' custom code that threw an error when one set of values were 0s; the
#' flowCore functions have better handling for this case.
#'
#' @param fcsFile A string containing the location of an .fcs file
#' @return Boolean: TRUE if there's an offset issue detected, FALSE otherwise
#'
#' @author Mario G Rosasco
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @export
checkForFcsByteOffsetIssue <- function(fcsFile){
  fileHandle <- file(fcsFile, "rb")
  # check for start/end validity using flowCore functions
  offsets <-
    tryCatch(
      flowCore:::findOffsets(fileHandle), 
      error = function(cond) {
        message(paste0("Byte offset issue in file ", fcsFile, ": ", cond, "\n"))
        return(NA)})
  close(fileHandle)
  
  if (length(offsets) == 1 && is.na(offsets)) {
    byteOffsetOk <- FALSE
  } else if (offsets[["datastart"]] > 0 && offsets[["dataend"]] > offsets[["datastart"]]) {
    byteOffsetOk <- TRUE
  } else stop(paste0("Error encountered checking for byte offset issue in file ", fcsFile, ".",
                     "Offset values are illogical.\n",
                     "Data start value: ", offsets[["datastart"]], "\n",
                     "Data end value: ", offsets[["dataend"]]))
  
  return(!byteOffsetOk)
}

#' Check the number of events in .fcs files
#'
#' Many clustering algorithms cannot handle data objects with no events, so this
#' utility is provided as a way to quickly check the number of events in a file.
#' The number of events is parsed from the text header using low-level file
#' I/O, so counting the events in very large .fcs files will not take
#' longer than smaller files. It was modified in Feb 2022 to use flowCore
#' functions that read the data location from the header tags if the values are
#' missing from the initial header string.
#'
#' @param fcsFile A string containing the location of an .fcs file
#' @return Numeric indicating the number of events listed in the .fcs header
#'
#' @author Mario G Rosasco
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @importFrom stringr str_extract str_replace str_replace_all regex
#' @export
getFcsNEvents <- function(fcsFile) {
  header <- flowCore::read.FCSheader(fcsFile)
  nEvents <- as.integer(header[[1]][["$TOT"]])
  return(nEvents)
}


#' Display the per-subject event counts in an experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment()}
#' @return A data frame containing the number of events detected for each subject
#'
#' @author Mario G Rosasco, Virginia Muir
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
#' @author Mario G Rosasco, Virginia Muir
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
#' @author Mario G Rosasco
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

#' Compute and return individual event intensities from a metaclustered experiment
#'
#' @param experiment A metaclustered discovrExperiment
#' @param cellSubsets A string or vector of strings indicating which cell subsets to return event z-scores from
#' @param metaclusters An integer or vector of integers indicating which metaclusters to return event z-scores from
#' @param subjects A string or vector of strings indicating which subjects to return event z-scores from
#' @param markers A string or vector of strings indicating which markers to return event z-scores from
#' @return A data frame containing arcsinh-transformed event intensities from the events satisfying the requested filters. For any filter field without a specified value, all possible events will be returned.
#'
#' @author Mario G Rosasco
#' @export
getEventIntensities <- function(experiment, cellSubsets = NA, metaclusters = NA, subjects = NA, markers = NA){
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

  message("Retrieving intensities for the selected events. Please be patient, this may take a moment...")

  # get metacluster indices
  metaclusterIndices <-
    data.frame(metacluster = experiment$colIndices) %>%
    tibble::rownames_to_column("sampRpClust") %>%
    dplyr::filter(.data$metacluster %in% metaclusters)

  # filter the data
  rowsToUse <- (experiment$mergedExpr$samp %in% subjects) & (experiment$mergedExpr$cellSubset %in% cellSubsets)
  colsToUse <- unique(c(markers, "samp", "cellSubset", "sampRpClust"))

  intensityData <-
    experiment$mergedExpr[rowsToUse, colsToUse] %>%
    # merge in metacluster indices and filter out metaclusters to exclude
    dplyr::inner_join(metaclusterIndices, by = "sampRpClust") %>%
    dplyr::rename(
      cluster = .data$sampRpClust,
      subject = .data$samp
    )

  return(intensityData)
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
#' @author Mario G Rosasco
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
    tidyr::gather("marker", "subjectStdDev", -.data$samp)

  message("Computing per-subject means...")
  parentMean <-
    experiment$mergedExpr %>%
    dplyr::select(-.data$cellSubset, -.data$RPclust, -.data$sampRpClust) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(.data$samp) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather("marker", "subjectMean", -.data$samp)

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
    # get filtered intensity data with metacluster info attached
    getEventIntensities(
      experiment = experiment,
      cellSubsets = cellSubsets,
      metaclusters = metaclusters,
      subjects = subjects,
      markers = markers
    ) %>%
    # track event ID to help re-spread table later
    dplyr::mutate(event = row_number()) %>%
    tidyr::gather("marker", "value", -.data$subject, -.data$cellSubset, -.data$cluster, -.data$event, -.data$metacluster) %>%
    # merge in stdev and mean
    dplyr::inner_join(parentMeanStdDev, by = c("subject" = "samp", "marker" = "marker")) %>%
    # compute z-score
    dplyr::mutate(zScore = (.data$value-.data$subjectMean)/.data$subjectStdDev) %>%
    # reformat data
    dplyr::select(.data$subject, .data$cellSubset, .data$marker, .data$zScore, .data$cluster, .data$metacluster, .data$event) %>%
    tidyr::spread("marker", "zScore") %>%
    dplyr::select(-.data$event)

  return(zScoreData)
}


#' Check a list of FCS files for a number of potential issues
#' 
#' Internal function used by \code{briDiscovr} to verify FCS file locations,
#' integrity, and optionally size relative to system memory. Split out as a
#' separate function because it is used by multiple other functions in this
#' package. Steps formerly contained with \code{setupDiscovrExperiment}.
#' 
#' @param fcsInfo A data frame containing file locations and information for
#' data to be included in a briDiscovr analysis. Typically loaded by another
#' \code{briDiscovr} function by specifying \code{fcsInfoFile}
#' @param checkMemory (default: TRUE) Logical indicating whether to check how
#' much system memory is available before loading the dataset. If TRUE, this
#' function will display a message and prevent data loading when the files take
#' up more than 80 percent of the available system memory.
#' @param verbose (default: TRUE) Logical specifying whether to display processing messages
#' @return logical, TRUE if the files in fcsInfo pass all checks. Function 
#' should exit with a \code{stop} statement if any checks are failed.
#' 
#' @seealso \code{\link{setupDiscovrExperiment}}, \code{\link{downsampleFcsList}}
#' @author Mario G Rosasco
#' @author Virginia Muir
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @export
#' @import flowCore
#' @import dplyr
#' @importFrom memuse Sys.meminfo mu.size
checkFcsFiles <- function(
    fcsInfo,
    checkMemory = TRUE,
    verbose = TRUE
){
  # Check for appropriate columns in fcsInfo
  if(!all(c("subject", "cellSubset", "filename") %in% names(fcsInfo)))
    stop("The file set as 'fcsInfoFile' must contain columns with names 'subject', 'cellSubset', and 'filename'.")
  
  # Check that fcs files can be found, and get size of all files
  filesThatDontExist <- c()
  sumFileSize <- 0
  for(currFile in fcsInfo$filename){
    if(!file.exists(currFile)){
      filesThatDontExist <- c(filesThatDontExist, currFile)
    } else {
      sumFileSize <- sumFileSize + file.info(currFile)$size
    }
  }
  if(length(filesThatDontExist > 0)){
    stop(
      "The following files could not be found:\n",
      paste0(filesThatDontExist, collapse = "\n")
    )
  }
  
  # Check that the files don't exceed the available system memory
  if(checkMemory){
    currMem <- memuse::Sys.meminfo()
    currSysMem <- mu.size(currMem$totalram, as.is = FALSE)
    if(sumFileSize > currSysMem){
      stop(
        "Total file size of dataset (", sumFileSize,
        " bytes) exceeds the total system memory (", currSysMem, " bytes)."
      )
    }
    currFreeMem <- mu.size(currMem$freeram, as.is = FALSE)
    if(checkMemory & (sumFileSize > 0.8*currFreeMem)){
      stop(
        "Total file size of dataset (", sumFileSize,
        " bytes) exceeds 80% of the available memory (", currFreeMem, " bytes). ",
        "To process the dataset anyway, you can set the 'checkMemory' option to FALSE, ",
        "but be aware this may cause issues on your system if you run out of memory."
      )
    }
    if(verbose){
      message("Total system memory: ", currSysMem, " bytes")
      message("Available system memory: ", currFreeMem, " bytes")
      message("Size of dataset on disk: ", sumFileSize, " bytes - OK")
    }
  }
  
  # check that each subject has only one .fcs file for each cellSubset
  duplicatedSamples <- duplicated(paste0(fcsInfo$subject, fcsInfo$cellSubset))
  if(any(duplicatedSamples)){
    stop(
      "The following subject/cell subset pairs were duplicated in the .fcs list:\n",
      apply(
        fcsInfo[duplicatedSamples, c("subject", "cellSubset")],
        1,
        function(r){paste(c(r[1], " - ", r[2], "\n"))}
      )
    )
  }
  
  # Check fcs files for byte offset issue that will prevent analysis
  filesWithOffsetIssues <- c()
  for (currFile in fcsInfo$filename) {
    if (checkForFcsByteOffsetIssue(currFile)) {
      filesWithOffsetIssues <- c(filesWithOffsetIssues,
                                 currFile)
    }
  }
  if (length(filesWithOffsetIssues > 0)) {
    stop(paste("Found", length(filesWithOffsetIssues), "files with byte offset issues. For details, see error message(s) above. Please try re-exporting these files:\n",
               paste0(filesWithOffsetIssues, collapse = "\n")))
  }
  
  # Check fcs files for 0 events, which will prevent analysis
  filesWithNoEvents <- c()
  for(currFile in fcsInfo$filename){
    if(!getFcsNEvents(currFile) > 0){
      filesWithNoEvents <- c(filesWithNoEvents, currFile)
    }
  }
  if(length(filesWithNoEvents > 0)){
    stop(paste(
      "The following files have no events and cannot be used. Please remove these from your list of .fcs files:\n",
      paste0(filesWithNoEvents, collapse = "\n")
    ))
  }
  
  return(TRUE)
}


#' Downsample a set of FCS files reproducibly, to make the event count more uniform and reduce memory usage
#'
#' This function downsamples a set of FCS cytometry data files to a maximum
#' event count. For each FCS file representing a population in 
#' \code{popsToDownsample}, it checks whether the file exceeds the maximum
#' event count, set by \code{maxEvents}. If the maximum event count is
#' exceeded, it selects a random set of rows. These rows can either be stored
#' as a list of vectors and passed to \code{setupDiscovrExperiment}, or new
#' FCS files can be written out to a specified directory, with a new fcsInfoFile
#' specifying the updated file paths. Results can be made reproducible by
#' passing a non-NULL value for \code{seed}.
#' @param fcsInfoFile a character string indicating the path to a file with
#' information on FCS files to be used.
#' @param popsToDownsample a character vector with the names of the cell
#' populations to be downsampled. Allows the user to keep all events from
#' certain populations even if the number of events exceed \code{maxEvents}, 
#' then \code{maxEvents} rows are selected from the FCS file for future use.
#' For example, users may want to keep all events from antigen-specific
#' populations, while downsampling the parent populations.
#' @param maxEvents numeric, the maximum number of events to be kept for each
#' FCS file. This cap is applied only to files matching populations in
#' \code{popsToDownsample}.
#' @param downsampleMode (default: "storeVectors") character, either 
#' "storeVectors" or "writeFiles" (or unique abbreviations thereof). For 
#' "storeVectors", a list of vectors is returned, with one element for each file
#' listed in \code{fcsInfoFile}; this should be passed to 
#' \code{setupDiscovrExperiment} as \code{downSampleRowList}. For "writeFiles", 
#' you must specify a valid directory in \code{outdirFcsFilesDownsampled} to 
#' write new FCS files to, and a valid path for \code{outfileFcsInfoDownsampled}.
#' @param outfileFcsInfoDownsampled character string, the location to write out 
#' the FCS info file with updated paths for the downsampled files.  Used only if
#' \code{downsampleMode == "writeFiles"}.
#' @param outdirFcsFilesDownsampled character string, the path for the directory
#' to write out the the downsampled files. Used only if 
#' \code{downsampleMode == "writeFiles"}.
#' @param seed (default: NULL) numeric, the seed to be passed to 
#' \code{set.seed} to make the downsampling reproducible.
#' @param verbose (default: TRUE) A logical specifying whether to display processing messages
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @export
#' @return Depends on \code{downsampleMode}. If 
#' \code{downsampleMode == "storeVectors"}, a list of vectors is returned, one 
#' for each row in the file at fcsInfoFile, with the row numbers of events to 
#' be included. This object should be passed to \code{setupDiscovrExperiment}
#' via argument \code{downsampleVectorList} to enforce the downsampling. If 
#' \code{downsampleMode == "writeFiles"}, \code{NULL} is returned.
#' 
#' @details Files specified in \code{fcsInfoFile} are checked for a number of
#' potential issues, including missing paths, no events, and byte offset issues.
downsampleFcsList <- function(
    fcsInfoFile, popsToDownsample = "all", maxEvents = 1E5,
    downsampleMode = "storeVectors",
    outfileFcsInfoDownsampled, outdirFcsFilesDownsampled,
    seed = NULL,
    verbose = TRUE
){
  # check and read in fcsInfo from file (adapted from setupDiscovrExperiment)
  if(!file.exists(fcsInfoFile))
    stop("The file set as 'fcsInfoFile' cannot be accessed. Please check the file name and path.")
  fcsInfo <- read.csv(fcsInfoFile, stringsAsFactors = FALSE)
  
  # run checks on this input and the FCS files
  checkFcsFiles(fcsInfo, checkMemory = FALSE, verbose = verbose)
  
  downsampleMode <- match.arg(downsampleMode, choices = c("storeVectors", "writeFiles"))
  # check that file output locations are valid
  if(downsampleMode == "writeFiles") {
    if(!dir.exists(dirname(outfileFcsInfoDownsampled)))
      stop("The output location specified in 'outfileFcsInfoDownsampled' is not valid. Please check that the containing directory exists and has write permissions.")
    if(!dir.exists(outdirFcsFilesDownsampled))
      stop("The output directory specified in 'outdirFcsFilesDownsampled' is not valid. Please check that this directory exists and has write permissions.")
  }
  
  # set popsToDownsample
  if(popsToDownsample == "all")
    popsToDownsample <- unique(fcsInfo$cellSubset)
  if(!all(popsToDownsample %in% fcsInfo$cellSubset))
    stop("Some populations specified in 'popsToDownsample' were not found in the FCS info file. Please verify that the populations to be downsampled match values in the 'cellSubset' of the FCS info file.")
  
  if(!is.null(seed)) set.seed(seed)
  
  # iterate over all files
  rowVectorList <- list()
  for (currRow in (1:nrow(fcsInfo))) {
    currFilePath <- fcsInfo$filename[currRow]
    # get number of events in file
    currFileNEvents <- briDiscovr::getFcsNEvents(currFilePath)
    # if all events to be kept, store vector of 1 to allEvents
    if((!fcsInfo$cellSubset[currRow] %in% popsToDownsample) | (currFileNEvents <= maxEvents)) {
      rowVectorList[[currRow]] <- 1:currFileNEvents
    } else {
      rowVectorList[[currRow]] <- 
        sample.int(n = currFileNEvents, size = maxEvents, replace = FALSE)
    }
  }
  
  # return/output results
  if(downsampleMode == "storeVectors") {
    return(rowVectorList)
  } else if (downsampleMode == "writeFiles") {
    # generate new FCS info file with updated path names
    fcsInfoDownsampled <- fcsInfo
    fcsInfoDownsampled$filename <-
      file.path(dirname(fcsInfoDownsampled$filename), paste0("downsampled_", basename(fcsInfoDownsampled$filename)))
    write.csv(fcsInfoDownsampled, file = outfileFcsInfoDownsampled, row.names = FALSE)
    
    # read in original FCS files, and write out downsampled FCS files
    for(currRow in (1:nrow(fcsInfo))) {
      currFcs <-
        flowCore::read.FCS(
          fcsInfo$filename[currRow],
          truncate_max_range = FALSE) # do not truncate max range for spectral flow
      if(!identical(rowVectorList[[currRow]], "all"))
        currFcs <- currFcs[rowVectorList[[currRow]],]
      write.FCS(currFcs, filename = fcsInfoDownsampled$filename[currRow])
    }
  }
  return()
}
