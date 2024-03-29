% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportedUtils.R
\name{checkFcsFiles}
\alias{checkFcsFiles}
\title{Check a list of FCS files for a number of potential issues}
\usage{
checkFcsFiles(fcsInfo, checkMemory = TRUE, verbose = TRUE)
}
\arguments{
\item{fcsInfo}{A data frame containing file locations and information for
data to be included in a briDiscovr analysis. Typically loaded by another
\code{briDiscovr} function by specifying \code{fcsInfoFile}}

\item{checkMemory}{(default: TRUE) Logical indicating whether to check how
much system memory is available before loading the dataset. If TRUE, this
function will display a message and prevent data loading when the files take
up more than 80 percent of the available system memory.}

\item{verbose}{(default: TRUE) Logical specifying whether to display processing messages}
}
\value{
logical, TRUE if the files in fcsInfo pass all checks. Function 
should exit with a \code{stop} statement if any checks are failed.
}
\description{
Internal function used by \code{briDiscovr} to verify FCS file locations,
integrity, and optionally size relative to system memory. Split out as a
separate function because it is used by multiple other functions in this
package. Steps formerly contained with \code{setupDiscovrExperiment}.
}
\seealso{
\code{\link{setupDiscovrExperiment}}, \code{\link{downsampleFcsList}}
}
\author{
Mario G Rosasco

Virginia Muir

Matthew J Dufort, \email{mdufort@benaroyaresearch.org}
}
