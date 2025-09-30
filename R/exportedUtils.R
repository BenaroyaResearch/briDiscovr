## Copyright (C) 2025  Matt Dufort, Mario Rosasco, Virginia Muir, and Benaroya Research Institute
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
  } else stop(paste0("Error encountered checking for byte offset issue in file ", fcsFile, ".\n",
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
    group_by(sample) %>%
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
    group_by(sample) %>%
    summarize(kClusters = max(RPclust)) %>%
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
#' @param cellSubsets A string or vector of strings indicating which cell subsets to return event intensities from
#' @param metaclusters An integer or vector of integers indicating which metaclusters to return event intensities from
#' @param subjects A string or vector of strings indicating which subjects to return event intensities from
#' @param markers A string or vector of strings indicating which markers to return event intensities from
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
    dplyr::filter(metacluster %in% metaclusters)

  # filter the data
  rowsToUse <- (experiment$mergedExpr$samp %in% subjects) & (experiment$mergedExpr$cellSubset %in% cellSubsets)
  colsToUse <- unique(c(markers, "samp", "cellSubset", "sampRpClust"))

  intensityData <-
    experiment$mergedExpr[rowsToUse, colsToUse] %>%
    # merge in metacluster indices and filter out metaclusters to exclude
    dplyr::inner_join(metaclusterIndices, by = "sampRpClust") %>%
    dplyr::rename(
      cluster = sampRpClust,
      subject = samp
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
    dplyr::select(-cellSubset, -RPclust, -sampRpClust) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(samp) %>%
    dplyr::summarise_all(sd) %>%
    tidyr::gather("marker", "subjectStdDev", -.data$samp)

  message("Computing per-subject means...")
  parentMean <-
    experiment$mergedExpr %>%
    dplyr::select(-cellSubset, -RPclust, -sampRpClust) %>%
    unique() %>% # remove duplicated rows; mergedExpr has both parent and gated
    dplyr::group_by(samp) %>%
    dplyr::summarise_all(mean) %>%
    tidyr::gather("marker", "subjectMean", -.data$samp)

  parentMeanStdDev <-
    inner_join(parentMean, parentStdDev, by = c("samp", "marker"))

  # get metacluster indices
  metaclusterIndices <-
    data.frame(metacluster = experiment$colIndices) %>%
    tibble::rownames_to_column("sampRpClust") %>%
    dplyr::filter(metacluster %in% metaclusters)

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
    dplyr::mutate(zScore = (value - subjectMean) / subjectStdDev) %>%
    # reformat data
    dplyr::select(subject, cellSubset, marker, zScore, cluster, metacluster, event) %>%
    tidyr::spread("marker", "zScore") %>%
    dplyr::select(-event)

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
  
  # Check fcs files for 0 events, which will prevent analysis
  # this step should be done prior to checking for the byte offset issue,
  # as files with 0 events will return an uninformative byte offset error message
  fcsInfo$nEvents <- sapply(fcsInfo$filename, getFcsNEvents)
  filesWithNoEvents <-
    fcsInfo$filename[fcsInfo$nEvents == 0]
  if(length(filesWithNoEvents > 0)){
    stop(paste(
      "The following files have no events and cannot be used. Please remove these from your list of .fcs files:\n",
      paste0(filesWithNoEvents, collapse = "\n")
    ))
  }
  
  # Check fcs files for subjects with too few events, which will prevent clustering
  # this is done at the "subject" level rather than the fcs file level because multiple cellSubsets are clustered together for the same subject
  # only send a warning for this, to allow for clustering approaches downstream that do not require a minimum number of events
  minEventsPerSubject <- 32 # 32 is the minimum number of events required for clustering with k=30
  subjectsWithTooFewEvents <-
    fcsInfo %>%
    group_by(subject) %>%
    summarize(nEvents = sum(nEvents)) %>%
    dplyr::filter(nEvents < minEventsPerSubject) %>%
    pull(subject)
  if(length(subjectsWithTooFewEvents > 0)){
    message(
      "The following subjects have fewer than ", minEventsPerSubject,
      " total cells across all cell subsets and cannot be clustered under the ",
      "default parameters for RPhenograph. Please remove these from your list ",
      "of .fcs files:\n",
      paste0(subjectsWithTooFewEvents, collapse = "\n")
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
#' @importFrom utils write.csv
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

#' Test a range of normalization options on a marker in a discovrExperiment
#' 
#' This function tests different normalization methods on a marker in a
#' \code{discovrExperiment} object. Its purpose is to facilitate selection of
#' the appropriate normalization method for markers where this is not
#' immediately apparent. It can be run on a set of markers and a set of 
#' different normalization methods to be tested on each marker. Alternatively,
#' it can be run with a single warpSet normalization method and a set of
#' different random number seeds, to evaluate the effect of the random number
#' seed on the warpSet normalization results. The results can be examined for
#' appropriate distributions across samples, or plotted sequentially using
#' \code{plotDensityNormalizedExprsDiscovrExperiment}
#' @param experiment A discovrExperiment created using
#' \code{clusterDiscovrExperiment}, \code{normalizeDiscovrExperiment}, or
#' \code{metaclusterDiscovrExperiment}. In order to return metacluster numbers
#' for each cell, its status must be "metaclustered".
#' @param markers a character vector indicating the markers on which to test
#' normalization
#' @param normalizationMethod (default: c("none", "zScore", "warpSet")) a
#' character vector with the normalization methods to be tested. All values
#' should be acceptable inputs for the argument 'normalizationMethod' in
#' \code{normalizeDiscovrExperiment}. If \code{seed} has length > 1, this
#' must be a single character string, either "warpSet" or "warpSet[#]".
#' @param seed (default: 12345) Numeric, the random number seed(s) for warpSet
#' normalization, passed to \code{normalizeWarpSetMergedExpr}. If testing
#' various random number seeds with warpSet normalization, this should be a
#' numeric vector of length > 1, and \code{normalizationMethod} should be
#' either "warpSet" or "warpSet[#]"
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @export
#' @return A reduced \code{discovrExperiment} object (or named list of objects),
#' each containing only the marker listed, with normalized values in
#' \code{experiment$mergedExprNormalizedScaled}. If \code{markers} has length >
#' 1, the function will return a named list, with each element containing a
#' \code{discovrExperiment} object normalized using the given method.
testNormalizationMethodByMarker <- function(
    experiment,
    markers,
    normalizationMethod = c("none", "zScore", "warpSet"),
    seed = 12345
){
  # check that the experiment is a discovrExperiment
  if(!is.discovrExperiment(experiment))
    stop("The input 'experiment' must be a discovrExperiment object.")
  
  # check that either normalizationMethod or seed is of length 1
  # and if length(seed) > 1, then normalizationMethod must be warpSet or warpSet[#]
  if(length(seed) > 1){
    if((length(normalizationMethod) > 1) | !all(str_detect(normalizationMethod, "^warpSet(?=[0-9]*)"))){
      stop("Multiple random number seeds can only run with a single warpSet normalization method")
    } else {
      message(
        paste0("Detected multiple random number seeds and normalization method ", normalizationMethod, ".\n",
               "Running ", normalizationMethod, " normalization with each random number seed."))
      testMode <- "multiSeed"
    }
  } else if(length(normalizationMethod) == 1){
    testMode <- "singleNormMethod"
  } else testMode <- "multiNormMethods"
  
  # reduce contents of input experiment object to contain only marker(s) of interest
  if(!is.character(markers)) stop("Input value for 'markers' must be a character vector")
  if(!all(markers %in% colnames(experiment$mergedExpr))) {
    stop(
      paste0("Some markers to normalize were not found in the provided experiment object:\n",
             paste(setdiff(markers, colnames(experiment$mergedExpr)), collapse=", ")))
  }
  experiment$markerInfo <-
    experiment$markerInfo[experiment$markerInfo$commonMarkerName %in% markers,]
  experiment$mergedExpr <-
    experiment$mergedExpr[
      , colnames(experiment$mergedExpr) %in% c(markers, "samp", "cellSubset", "RPclust", "sampRpClust")]
  experiment$clusteringMarkers <- markers
  experiment$markersNormalized <- markers
  # remove existing normalization information if present
  experiment$mergedExprNormalizedScaled <- NULL
  experiment$clusterMeansNormalizedScaled <- NULL
  
  # run the normalization based on modes determined above
  if(testMode == "singleNormMethod"){
    experimentNormalized <-
      normalizeDiscovrExperiment(
        experiment,
        normalizationInfo = normalizationMethod,
        seed = seed,
        verbose = FALSE)
  } else if(testMode == "multiNormMethods"){
    experimentNormalized <- list()
    for(normMethod.tmp in normalizationMethod){
      experimentNormalized[[normMethod.tmp]] <-
        normalizeDiscovrExperiment(
          experiment,
          normalizationInfo = normMethod.tmp,
          seed = seed,
          verbose = FALSE)
    }
  } else if(testMode == "multiSeed"){
    experimentNormalized <- list()
    for(seed.tmp in seed){
      experimentNormalized[[paste0(normalizationMethod, "_seed", seed.tmp)]] <-
        normalizeDiscovrExperiment(
          experiment,
          normalizationInfo = normalizationMethod,
          seed = seed.tmp,
          verbose = FALSE)
    }
  }
  return(experimentNormalized)
}

#' Calculate UMAP coordinates for the cells from a discovrExperiment object 
#'
#' This function generates a UMAP from the cells in a discovrExperiment object.
#' To do this, either uses normalized data in the discovrExperiment object or
#' applies normalization to each marker similar to the briDiscovr normalization
#' / metaclustering process, optionally downsamples the cells to speed the
#' process (with downsampling frequency tunable at the cell population level),
#' and then runs the UMAP algorithm. The UMAP algorithm is run using the
#' \code{umap} package, which is a wrapper for the \code{uwot} package. Results
#' can be made reproducible by passing a non-NULL value for \code{seed}. This
#' function returns a data frame with the UMAP coordinates for each cell, as
#' well as the original cell population, sample information, and metacluster if
#' available. The outputs are intended to be visualized using plotting software
#' such as ggplot2.
#' @param experiment A discovrExperiment created using
#' \code{setupDiscovrExperiment}, \code{clusterDiscovrExperiment},
#' \code{normalizeDiscovrExperiment}, or \code{metaclusterDiscovrExperiment}.
#' In order to return metacluster numbers for each cell, its status must be
#' "metaclustered". If normalization has not already been run, it will be run
#' using \code{normalizeDiscovrExperiment} prior to running UMAP.
#' @param umapMarkers A character vector, the markers to be used for UMAP. For
#' the default value, NULL, the function extracts the set of markers from the
#' "clusteringMarkers" element of the discovrExperiment object.
#' @param downsampleBy character, specifying which downsampling method to use.
#' Acceptable values are "frequency", "number", or unique partial matches.
#' Defaults to "frequency" to replicate behavior of earlier versions. With
#' method "frequency", the cells from each sample are trimmed to keep 1 cell
#' for every N cells in the sample. With method "number", the cells are trimmed
#' to keep up to N cells from each sample. If "frequency" is used, the
#' "downsampleFreq" parameter is required and "downsampleNumber" will be
#' ignored. Similarly, if "number" is used, the "downsampleNumber" parameter is
#' required and "downsampleFreq" will be ignored.
#' @param downsampleFreq numeric, specifying how to downsample the cells prior
#' to running UMAP. This approach keeps 1 cell for every N total cells in each
#' population in each sample. Several alternative methods can be used by providing
#' different numeric vectors. If a single value is provided, all
#' populations are downsampled to this frequency. If a vector of length 2 is
#' provided (optionally with elements named "parentPopulation" and
#' "childPopulations"), the "parentPopulation" or first element is used as the
#' frequency for the parent population (extracted from the discovrExperiment
#' object), and the "childPopulations" or second element is used as the
#' frequency for the child populations. If a named vector is provided, the names
#' must match the cell populations, and the values are the frequencies
#' to downsample each population to. If NULL, no downsampling is performed. The
#' default is c("parentPopulation" = 100, "childPopulations" = 1), which retains
#' all cells from child populations and subsets the parent population to 1/100.
#' Note that downsampling is based on the order of the cells in
#' discovrExperiment, so changes that alter the order of cells will make the
#' downsampling results non-reproducible.
#' @param downsampleNumber numeric, specifying how to downsample the cells
#' prior to running UMAP. This approach keeps up to N cells for each population
#' in each sample. Several alternative methods can be used by providing
#' different numeric vectors. If a single value is provided, all populations
#' are downsampled to this number of cells. If a vector of length 2 is provided
#' (optionally with elements named "parentPopulation" and "childPopulations"),
#' the "parentPopulation" or first element is used as the number for the parent
#' population (extracted from the discovrExperiment object), and the
#' "childPopulations" or second element is used as the number for the child
#' populations. If a named vector is provided, the names must match the cell
#' populations, and the values are the numbers to downsample each population
#' to. If NULL, no downsampling is performed. To keep all cells in a population,
#' set to Inf. The default is c("parentPopulation" = 100, "childPopulations" =
#' Inf), which retains all cells from child populations and subsets the parent
#' population to 100 cells. Note that downsampling is based on the order of the
#' cells in discovrExperiment, so changes that alter the order of cells will
#' make the downsampling results non-reproducible.
#' @param normalizationInfo (default: NULL) character string or vector
#' specifying normalization method to be applied, passed to
#' \code{normalizeDiscovrExperiment}. NOTE: if this argument is non-null,
#' the specified normalization will be applied to the input 'experiment', even
#' if that object has already had normalization applied. Thus, including a
#' non-null value here will override any existing normalization prior to
#' running UMAP.
#' @param seed (default: NULL) numeric, the seed to be passed to 
#' \code{set.seed} to make the UMAP (more) reproducible. If NULL, no seed is set.
#' @param returnUmapObject (default: FALSE) logical, if TRUE, returns the full
#' UMAP output object and the data frame of cell information as a list. If
#' FALSE, returns only the data frame of cell information.
#' @param returnExpressionNormalizedScaled (default: FALSE) logical, if TRUE,
#' includes the normalized, scaled expression values for the markers used in the
#' UMAP in the output data frame. If FALSE, the data frame only contains the
#' UMAP coordinates and cell information.
#' @param ... optional arguments passed to \code{umap::umap}.
#' @importFrom umap umap
#' @import dplyr
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @export
#' @return A data frame containing the UMAP coordinates for each cell, as columns
#' 'UMAP1' and 'UMAP2', and the original cell population and sample information.
#' The data frame also contains the metacluster information, if available.
#' If \code{returnUmapObject} is TRUE, returns a list with the data frame of
#' cell information as element 'data' and the UMAP output object as element
#' 'umapObject'.
#' If \code{returnExpressionNormalizedScaled} is TRUE, the data frame also
#' contains the normalized expression values for the markers used in the UMAP.
runUmapDiscovrExperiment <- function(
    experiment,
    umapMarkers = NULL,
    downsampleBy = "frequency",
    downsampleFreq = c("parentPopulation" = 100, "childPopulations" = 1),
    downsampleNumber = c("parentPopulation" = 100, "childPopulations" = Inf),
    normalizationInfo = NULL,
    seed = NULL,
    returnUmapObject = FALSE,
    returnExpressionNormalizedScaled = FALSE,
    ...
){
  # check that the experiment is a discovrExperiment
  if(!is.discovrExperiment(experiment))
    stop("The input 'experiment' must be a discovrExperiment object.")
  
  # normalize expression values if normalizationInfo is non-null, even if normalization has already been done
  if(!is.null(normalizationInfo)) {
    message("Detected normalization specifications in 'normalizationInfo'. Running normalization prior to UMAP.\n")
    experiment <- normalizeDiscovrExperiment(experiment, normalizationInfo = normalizationInfo, ...)
    # check that normalized expression values are present, or normalizationInfo is provided
  } else if(experiment$status %in% c("normalized", "metaclustered")) {
    if(is.null(experiment$mergedExprNormalizedScaled))
      stop("The input 'experiment' has status", experiment$status,
           " but does not include normalized expression values. Please ensure that it is run through 'normalizeDiscovrExperiment'")
  } else if(experiment$status %in% c("initialized", "clustered")) {
    if(is.null(experiment$markerInfo$normalizationMethod)) {
      stop("The input 'experiment' has status ", experiment$status, ".\n",
           "In order to run this function on a 'clustered' discovrExperiment object,\n",
           "the input 'experiment' must include normalization specifications in 'markerInfo'\n",
           "or you must provide a non-null input for 'normalizationInfo' that can be passed to 'normalizedDiscovrExperiment'.")
    } else {
      message(
        "The input 'experiment' has not yet had marker expression normalized using the\n",
        "'normalizeDiscovrExperiment' function. Normalization will now be performed prior to metaclustering,\n",
        "using the normalizationMethod specifications provided in 'markerInfo'.\n",
        "Note that this function does not perform default normalization unless explicitly told to do so\n",
        "using the argument 'normalizationInfo', in order to avoid unintended results.")
      
      experiment$status <- "clustered" # set status to clustered to avoid some hickups with normalizeDiscovrExperiment
      experiment <- normalizeDiscovrExperiment(experiment)
    }
  }

# check that umapMarkers is valid, and set it to the default if NULL
  if(is.null(umapMarkers)) {
    umapMarkers <- experiment$clusteringMarkers
  } else if(!all(umapMarkers %in% colnames(experiment$mergedExprNormalizedScaled))) {
    stop(
      paste0(
        "The input 'umapMarkers' specifies markers that are not present in the expression data\n",
        "('mergedExpr' or 'mergedExprNormalizedScaled') of the discovrExperiment object.\n",
        "Please verify the inclusion of the following markers: ",
        paste(setdiff(umapMarkers, colnames(experiment$mergedExprNormalizedScaled)), 
              collapse = ", ")))
  } else if(length(umapMarkers) < 2) {
    stop("The input 'umapMarkers' must contain at least two markers.")
  } else umapMarkers <- unique(umapMarkers)
  
  # collapse cellSubset field to a character vector (rather than lists)
  if (!any(lengths(experiment$cellSubset) > 1)) {
    # if all contain only a single element, collapse to a character vector
    experiment$mergedExprNormalizedScaled$cellSubset <-
      sapply(experiment$mergedExprNormalizedScaled$cellSubset, `[[`, 1)
  } else {
    # if any contain multiple elements, use the first non-parent population
    experiment$mergedExprNormalizedScaled$cellSubset <-
      lapply(
        experiment$mergedExprNormalizedScaled$cellSubset,
        \(x) {
          if (length(x) == 1) x else x[which(x != experiment$parentPopulation)[1]]}) %>%
      unlist()
  }

  # check that downsampleBy matches expected values
  downsampleBy <- match.arg(downsampleBy, choices = c("frequency", "number"))

  # fork downsampling based on downsampleBy
  if (downsampleBy == "frequency") {
    # run downsampling using downsampleFreq
    
    # check that downsampleFreq contains integers > 0
    if(!all(downsampleFreq %% 1 == 0) | any(downsampleFreq <= 0))
      stop("The input 'downsampleFreq' must consist of positive integers")
    
    # check that downsampleFreq is valid, and structure it as a named vector
    if(is.null(downsampleFreq)) {
      # if downsampleFreq is NULL, no downsampling is performed
      downsampleFreq <-
        c("parentPopulation" = 1, "childPopulations" = 1)
    } else if(!is.null(names(downsampleFreq))) {
      if (!(setequal(names(downsampleFreq), c("parentPopulation", "childPopulations")) |
            setequal(names(downsampleFreq), experiment$mergedExprNormalizedScaled$cellSubset)))
        stop(
          paste0(
            "If the input 'downsampleFreq' is a named vector, the names must ",
            "match the values of cellSubset in the discovrExperiment object, or ",
            "be named 'parentPopulation' and 'childPopulations'."))
    } else if(is.numeric(downsampleFreq) && length(downsampleFreq) == 1) {
      # if downsampleFreq is of length 1, it is used for both parent and child
      downsampleFreq <-
        c("parentPopulation" = downsampleFreq, "childPopulations" = downsampleFreq)
    } else if(is.numeric(downsampleFreq) && length(downsampleFreq) == 2) {
      # if downsampleFreq is of length 2, it is used for parent and child
      if(setequal(names(downsampleFreq), c("parentPopulation", "childPopulations")))
        downsampleFreq <-
          c("parentPopulation" = downsampleFreq[["parentPopulation"]],
            "childPopulations" = downsampleFreq[["childPopulations"]])
      downsampleFreq <-
        c("parentPopulation" = downsampleFreq[[1]], "childPopulations" = downsampleFreq[[2]])
    } else if(is.numeric(downsampleFreq) && length(downsampleFreq) > 2) {
      # if downsampleFreq is of length > 2, it is used for each population
      if(length(downsampleFreq) != length(unique(names(downsampleFreq))))
        stop("The input 'downsampleFreq' contains duplicate names. Please provide a named vector with unique names.")
      if(!setequal(names(downsampleFreq), experiment$mergedExprNormalizedScaled$cellSubset))
        stop(
          "The input 'downsampleFreq' does not contain the same names as the cell populations in the\n",
          "'mergedExpr' or 'mergedExprNormalizedScaled' elements of the discovrExperiment object.")
    } else stop("The input 'downsampleFreq' is not valid. Please provide NULL or a numeric vector of length > 0.")
    
    # standardize downsampleFreq to be a named vector with the same names as the cell populations
    if(setequal(names(downsampleFreq), c("parentPopulation", "childPopulations"))) {
      downsampleFreq <-
        c(downsampleFreq[["parentPopulation"]],
          rep(downsampleFreq[["childPopulations"]],
              length(unique(experiment$mergedExprNormalizedScaled$cellSubset)) - 1))
      names(downsampleFreq) <-
        c(experiment$parentPopulation,
          setdiff(experiment$mergedExprNormalizedScaled$cellSubset, experiment$parentPopulation))
    }
  
    
  } else if(downsampleBy == "number") {
    # run downampling by downsampleNumber
    
    # check that downsampleNumber contains integers > 0 or Inf
    if(!all((downsampleNumber %% 1 == 0) | (downsampleNumber == Inf)) | any(downsampleNumber <= 0))
      stop("The input 'downsampleNumber' must consist of positive integers or infinity")
    
    # check that downsampleNumber is valid, and structure it as a named vector
    if(is.null(downsampleNumber)) {
      # if downsampleNumber is NULL, no downsampling is performed
      downsampleNumber <-
        c("parentPopulation" = Inf, "childPopulations" = Inf)
    } else if(!is.null(names(downsampleNumber))) {
      if (!(setequal(names(downsampleNumber), c("parentPopulation", "childPopulations")) |
            setequal(names(downsampleNumber), experiment$mergedExprNormalizedScaled$cellSubset)))
        stop(
          paste0(
            "If the input 'downsampleNumber' is a named vector, the names must ",
            "match the values of cellSubset in the discovrExperiment object, or ",
            "be named 'parentPopulation' and 'childPopulations'."))
    } else if(is.numeric(downsampleNumber) && length(downsampleNumber) == 1) {
      # if downsampleNumber is of length 1, it is used for both parent and child
      downsampleNumber <-
        c("parentPopulation" = downsampleNumber, "childPopulations" = downsampleNumber)
    } else if(is.numeric(downsampleNumber) && length(downsampleNumber) == 2) {
      # if downsampleNumber is of length 2, it is used for parent and child
      if(setequal(names(downsampleNumber), c("parentPopulation", "childPopulations")))
        downsampleNumber <-
          c("parentPopulation" = downsampleNumber[["parentPopulation"]],
            "childPopulations" = downsampleNumber[["childPopulations"]])
      downsampleNumber <-
        c("parentPopulation" = downsampleNumber[[1]], "childPopulations" = downsampleNumber[[2]])
    } else if(is.numeric(downsampleNumber) && length(downsampleNumber) > 2) {
      # if downsampleNumber is of length > 2, it is used for each population
      if(length(downsampleNumber) != length(unique(names(downsampleNumber))))
        stop("The input 'downsampleNumber' contains duplicate names. Please provide a named vector with unique names.")
      if(!setequal(names(downsampleNumber), experiment$mergedExprNormalizedScaled$cellSubset))
        stop(
          "The input 'downsampleNumber' does not contain the same names as the cell populations in the\n",
          "'mergedExpr' or 'mergedExprNormalizedScaled' elements of the discovrExperiment object.")
    } else stop("The input 'downsampleNumber' is not valid. Please provide NULL or a numeric vector of length > 0.")
    
    # standardize downsampleNumber to be a named vector with the same names as the cell populations
    if(setequal(names(downsampleNumber), c("parentPopulation", "childPopulations"))) {
      downsampleNumber <-
        c(downsampleNumber[["parentPopulation"]],
          rep(downsampleNumber[["childPopulations"]],
              length(unique(experiment$mergedExprNormalizedScaled$cellSubset)) - 1))
      names(downsampleNumber) <-
        c(experiment$parentPopulation,
          setdiff(experiment$mergedExprNormalizedScaled$cellSubset, experiment$parentPopulation))
    }
  }
  
  # extract the cell populations and the merged, normalized, scaled expression matrix
  # this is run after the initial downsampling steps so that it doesn't run if there's anything wrong with the downsampling inputs
  exprData <-
    experiment$mergedExprNormalizedScaled %>%
    dplyr::select(samp, cellSubset, any_of("sampRpClust"), all_of(umapMarkers))
  
  # downsample the data
  
  if(downsampleBy == "frequency") {
    # run downampling by downsampleFreq
  
    # downsample the data, keeping 1 cell per downsamplefreq
    if(any(downsampleFreq > 1)) {
      exprData <-
        exprData %>%
        left_join(
          data.frame(
            cellSubset = names(downsampleFreq),
            freq = downsampleFreq),
          by = "cellSubset") %>%
        group_by(samp, cellSubset) %>%
        # keep 1 cell per freq; use remainder of 1 so that the first cell is retained if freq > 1
        mutate(sampCellSubsetIter = 1:n(),
               freqRemainder = sampCellSubsetIter %% freq,
               keepCell = case_when(freq == 1 ~ TRUE,
                                    freqRemainder == 1 ~ TRUE,
                                    freq > 1 & freqRemainder != 1 ~ FALSE)) %>% 
        dplyr::filter(keepCell == TRUE) %>%
        ungroup() %>%
        select(-freq, -freqRemainder, -sampCellSubsetIter, -keepCell)
    }
    
  } else if(downsampleBy == "number") {
    # run downampling by downsampleNumber
    
    # downsample the data, keeping N cells per population/sample
    if(any(downsampleNumber < Inf)) {
      # create a data frame with the rows to keep (a little clunky but much more efficient than other options I tried)
      exprDataCellsToKeepBySampCellSubset <-
        exprData %>%
        left_join(
          data.frame(
            cellSubset = names(downsampleNumber),
            nCellsToKeep = downsampleNumber),
          by = "cellSubset") %>%
        group_by(samp, cellSubset, nCellsToKeep) %>%
        # set nCellsToKeep to be not more than the number of cells
        summarize(nCellsBySampCellSubset = n()) %>%
        mutate(nCellsToKeep = min(nCellsToKeep, nCellsBySampCellSubset)) %>%
        ungroup()
      exprDataCellsToKeepBySampCellSubset[["rowNumberBySampCellSubset"]] <-
        mapply(\(nCellsBySampCellSubset, nCellsToKeep) round(seq(from = 1, to = nCellsBySampCellSubset, length.out = nCellsToKeep)),
               exprDataCellsToKeepBySampCellSubset$nCellsBySampCellSubset, exprDataCellsToKeepBySampCellSubset$nCellsToKeep)
      exprDataCellsToKeepBySampCellSubset <-
        unnest(exprDataCellsToKeepBySampCellSubset, cols = rowNumberBySampCellSubset) %>%
        mutate(keepCell = TRUE)

      # filter the expression data
      exprData <-
        exprData %>%
        group_by(samp, cellSubset) %>%
        mutate(rowNumberBySampCellSubset = row_number()) %>%
        ungroup() %>%
        # merge in rows to keep
        left_join(exprDataCellsToKeepBySampCellSubset, by = c("samp", "cellSubset", "rowNumberBySampCellSubset")) %>%
        dplyr::filter(keepCell == TRUE) %>%
        select(-rowNumberBySampCellSubset, -nCellsToKeep, -nCellsBySampCellSubset, -keepCell)
    }
  }
  
  # run UMAP on the downsampled data
  if(!is.null(seed)) set.seed(seed)
  umapRes <-
    umap::umap(
      exprData %>%
        select(all_of(umapMarkers)) %>%
        as.matrix(),
      ...)
  
  # if discovrExperiment object has been metaclustered, add metaclusters to exprData
  if(experiment$status == "metaclustered") {
    exprData <-
      exprData %>%
      left_join(
        data.frame(
          sampRpClust = names(experiment$colIndices),
          metacluster = unname(experiment$colIndices)),
        by = "sampRpClust") %>%
      select(-sampRpClust)
  }
  
  # remove expression values from the data frame, unless returnExpressionNormalizedScaled is TRUE
  if(!isTRUE(returnExpressionNormalizedScaled))
    exprData <- exprData %>%
    select(samp, cellSubset, any_of("metacluster"))
  
  # add the UMAP coordinates to the data frame
  exprData <-
    exprData %>%
    bind_cols(
      as.data.frame(umapRes$layout) %>%
        setNames(paste0("UMAP", seq_len(ncol(umapRes$layout))))) %>%
    as.data.frame()
  
  if(isTRUE(returnUmapObject)) {
    # return the UMAP object
    return(
      list(data = exprData,
           umapObject = umapRes))
  } else {
    # return the data frame with UMAP coordinates
    return(exprData)
  }
}

#' Plot normalized vs. pre-normalization marker level distributions for a discovrExperiment object
#'
#' This function generates a set of density plots of the pre- and
#' post-normalization expression values from a discovrExperiment object. If
#' the object has status "normalized" or "metaclustered", the normalized values
#' will be used. If the object has status "clustered", normalization will be
#' attempted using information provided, passed to
#' \code{normalizeDiscovrExperiment}. The function returns a list of plots as
#' ggplot2 objects, and can optionally output plots to a single file containing
#' a separate page for each marker. The purpose of this function is to evaluate
#' potential normalization methods for each marker, based on distributions of
#' the marker levels before and after normalization.
#' @param experiment A discovrExperiment created using
#' \code{clusterDiscovrExperiment}, \code{normalizeDiscovrExperiment},
#' \code{metaclusterDiscovrExperiment}. Must have status "clustered",
#' "normalized", or "metaclustered".
#' @param clusteringMarkersOnly (default: NULL)
#' @param normalizationInfo (default: NULL) Optional object containing
#' information on how to normalize, passed to \code{normalizeDiscovrExperiment}
#' (see documentation of that function for usage details). If this argument
#' is non-null, normalization will be performed prior to plotting, even if the
#' 'experiment' object has already been normalized. This allows for testing a
#' range of normalization approaches on a discovrExperiment object. A non-null
#' value must be provided if \code{experiment} has status "clustered" or does
#' not contain the element \code{mergedExprNormalizedScaled} which is generated
#' by \code{normalizeDiscovrExperiment}.
#' @param filenameOut (default: NULL) Path for outputing a file with plots. If
#' null, no file is generated. The extension is used to determine plot type,
#' which can be "pdf", "png", or "svg".
#' @param width (default: 12) Numeric, width of output plot device in inches
#' @param height (default: 6) Numeric, height of output plot device in inches
#' @param ... optional arguments passed to \code{normalizeDiscovrExperiment}
#' @import ggplot2
#' @importFrom stringr str_extract str_to_lower
#' @importFrom rlang sym :=
#' @importFrom grDevices jpeg pdf svg tiff
#' @author Matthew J Dufort, \email{mdufort@@benaroyaresearch.org}
#' @export
#' @return A list of ggplot objects, one for each marker. Plots will be output
#' to file if filenameOut is specified. To suppress output to the default
#' plotting device, assign the result to an object, or run the function with
#' \code{invisible()}
plotDensityNormalizedExprsDiscovrExperiment <- function(
    experiment,
    clusteringMarkersOnly = TRUE,
    normalizationInfo = NULL,
    filenameOut = NULL,
    width = 12, height = 6,
    ...
){
  # check that the experiment is a discovrExperiment
  if(!is.discovrExperiment(experiment))
    stop("The input 'experiment' must be a discovrExperiment object.")
  
  # check status of the experiment
  if(!(experiment$status %in% c("clustered", "normalized", "metaclustered")))
    stop("The input 'experiment' must have status 'clustered', 'normalized', 'metaclustered' in order to use this function. ",
         "The provided 'experiment' object has status ", experiment$status)
  
  # # normalize expression values if normalizationInfo is non-null
  # if(!is.null(normalizationInfo)) {
  #   experiment <- normalizeDiscovrExperiment(experiment, normalizationInfo = normalizationInfo, ...)
  # # check that normalized expression values are present, or normalizationInfo is provided
  # } else if(experiment$status %in% c("normalized", "metaclustered")) {
  #   if(is.null(experiment$mergedExprNormalizedScaled))
  #     stop("The input 'experiment' has status", experiment$status, " but does not include normalized expression values.\n",
  #          "Please ensure that it is run through 'normalizeDiscovrExperiment'")
  # } else if(experiment$status == "clustered") {
  #   stop("The input 'experiment' has status 'clustered'. In order to run this function on a 'clustered' discovrExperiment object, ",
  #        "you must provide a non-null input for 'normalizationInfo' that can be passed to 'normalizedDiscovrExperiment'.")
  # }
  
  # determine markers to plot
  if(isTRUE(clusteringMarkersOnly)) {
    markersToPlot <- experiment$clusteringMarkers
  } else {
    markersToPlot <-
      setdiff(
        intersect(colnames(experiment$mergedExpr), colnames(experiment$mergedExprNormalizedScaled)),
        c("samp", "cellSubset", "RPclust"))
  }
  
  # generate plots
  plotList <- list()
  for (marker.tmp in markersToPlot) {
    plotList[[marker.tmp]] <-
      ggplot(
        bind_rows(
          experiment$mergedExpr %>%
            dplyr::select(samp, markerExpr := !!rlang::sym(marker.tmp)) %>%
            mutate(normStatus = "arcsinh-transformed"),
          experiment$mergedExprNormalizedScaled %>%
            dplyr::select(samp, markerExpr := !!rlang::sym(marker.tmp)) %>%
            mutate(normStatus = "normalized")),
        mapping = aes(x = markerExpr, group = samp)) +
      geom_density() +
      facet_wrap(~normStatus, scales = "free") +
      labs(
        x = "marker level", y = "density",
        title = 
          paste0(
            marker.tmp,
            " (normalization: ",
            experiment$markerInfo$normalizationMethod[match(marker.tmp, experiment$markerInfo$commonMarkerName)],
            ")"))
  }
  

  # output plots if specified to do so
  if(!is.null(filenameOut)) {
    # check extension and open plotting device
    extensionOut <- str_extract(str_to_lower(filenameOut), "(?<=\\.)(pdf|png|jpeg|tiff|svg)(?=$)")
    if(is.na(extensionOut)) stop("Filename specified for output must have extension 'pdf', 'png', 'jpeg', 'tiff', or 'svg'")
        if(extensionOut == "pdf") {
      pdf(file = filenameOut, width = width, height = height)
    } else if(extensionOut == "png") {
      png(filename = filenameOut, width = width, height = height, units = "inches")
    } else if(extensionOut == "jpeg") {
      jpeg(filename = filenameOut, width = width, height = height, units = "inches")
    } else if(extensionOut == "tiff") {
      tiff(filename = filenameOut, width = width, height = height, units = "inches")
    } else if(extensionOut == "svg") {
      svg(filename = filenameOut, width = width, height = height)
    }
    
    for (plot.tmp in plotList) print(plot.tmp)
    
    invisible(dev.off())
  }
  
  return(plotList)
}
