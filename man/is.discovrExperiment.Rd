% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discovrExperiment.R
\name{is.discovrExperiment}
\alias{is.discovrExperiment}
\alias{discovrExperiment}
\title{discovrExperiment S3 class}
\usage{
## S3 methods for class 'discovrExperiment'
is.discovrExperiment(x)
}
\arguments{
\item{x}{object to be tested for class discovrExperiment}
}
\description{
Functions to create and manage discovrExperiment objects.
}
\details{
The discovrExperiment is implemented as a handy way to pass a lot of DISCOV-R data
around between the setup, clustering, and metaclustering steps. In normal usage,
users should not need to interact directly with the underlying data structure, and
should rely on the setup, clustering, and metaclustering functions implemented in the
\code{briDiscovr} package to manage discovrExperiment objects.
}
\seealso{
\code{\link{setupDiscovrExperiment}}
}
\author{
Mario G Rosasco
}
