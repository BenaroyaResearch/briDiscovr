## Copyright (C) 2020  Mario Rosasco, Virginia Muir, and Benaroya Research Institute
##
## This file is part of the briDiscovr package

#' discovrExperiment S3 class
#'
#' Functions to create and manage discovrExperiment objects.
#'
#' The discovrExperiment is implemented as a handy way to pass a lot of DISCOV-R data
#' around between the setup, clustering, and metaclustering steps. In normal usage,
#' users should not need to interact directly with the underlying data structure, and
#' should rely on the setup, clustering, and metaclustering functions implemented in the
#' \code{briDiscovr} package to manage discovrExperiment objects.
#'
#' @usage
#'
#' ## S3 methods for class 'tcrGraph'
#' is.discovrExperiment(x)
#'
#' @param x object to be tested for class discovrExperiment
#'
#' @aliases discovrExperiment is.discovrExperiment
#'
#' @seealso \code{\link{setupDiscovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
is.discovrExperiment <- function(x){ inherits(x, "discovrExperiment") }

#' Print method for discovrExperiment objects
#'
#' Prints summary information about the experiment structure.
#'
#' @param x discovrExperiment object
#' @param ... other arguments to be passed to generic print method
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @method print discovrExperiment
#' @export
print.discovrExperiment <- function(x, ...){
  if(!is.discovrExperiment(x))
    stop("The argument to print.discovrExperiment must be an object of type 'discovrExperiment'.")
  cat(
    sprintf(
      "An object of class 'discovrExperiment'\nAll markers: %s\nClustering markers: %s\nExperiment status: %s\n",
      paste0(x$markerInfo$commonMarkerName, collapse = ", "),
      paste0(x$clusteringMarkers, collapse = ", "),
      x$status
    )
  )
  if(length(list(...)) > 0){
    warning("Additional arguments were passed to print.discovrExperiment and were ignored.")
  }
}

#' Get all sample IDs in an experiment
#'
#' @param x discovrExperiment object
#'
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
getSampleNames <- function(x){
  if(!is.discovrExperiment(x))
    stop("The argument to 'getSampleNames' must be an object of type 'discovrExperiment'.")
  return(unique(x$mergedExpr$samp))
}
