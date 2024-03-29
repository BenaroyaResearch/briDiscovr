% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportedUtils.R
\name{getEventIntensities}
\alias{getEventIntensities}
\title{Compute and return individual event intensities from a metaclustered experiment}
\usage{
getEventIntensities(
  experiment,
  cellSubsets = NA,
  metaclusters = NA,
  subjects = NA,
  markers = NA
)
}
\arguments{
\item{experiment}{A metaclustered discovrExperiment}

\item{cellSubsets}{A string or vector of strings indicating which cell subsets to return event z-scores from}

\item{metaclusters}{An integer or vector of integers indicating which metaclusters to return event z-scores from}

\item{subjects}{A string or vector of strings indicating which subjects to return event z-scores from}

\item{markers}{A string or vector of strings indicating which markers to return event z-scores from}
}
\value{
A data frame containing arcsinh-transformed event intensities from the events satisfying the requested filters. For any filter field without a specified value, all possible events will be returned.
}
\description{
Compute and return individual event intensities from a metaclustered experiment
}
\author{
Mario G Rosasco
}
