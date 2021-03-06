% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportedUtils.R
\name{getFcsNEvents}
\alias{getFcsNEvents}
\title{Check the number of events in .fcs files}
\usage{
getFcsNEvents(fcsFile)
}
\arguments{
\item{fcsFile}{A string containing the location of an .fcs file}
}
\value{
Numeric indicating the number of events listed in the .fcs header
}
\description{
Many clustering algorithms cannot handle data objects with no events, so this
utility is provided as a way to quickly check the number of events in a file.
The number of events is parsed from the text header using low-level file
I/O, so counting the events in very large .fcs files will not take
longer than smaller files.
}
\author{
Mario G Rosasco, \email{mrosasco@benaroyaresearch.org}
}
