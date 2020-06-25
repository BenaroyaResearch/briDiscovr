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
#' @aliases is.discovrExperiment
#'
#' @seealso \code{\link{setupDiscovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
is.discovrExperiment <- function(x){ inherits(x, "discovrExperiment") }
