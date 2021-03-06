% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{jaccard_coeff}
\alias{jaccard_coeff}
\title{Compute jaccard coefficient between nearest-neighbor sets}
\usage{
jaccard_coeff(idx)
}
\arguments{
\item{idx}{NumericMatrix}
}
\value{
NumericMatrix
}
\description{
Weights of both i->j and j->i are recorded if they have intersection. In this case
w(i->j) should be equal to w(j->i). In some case i->j has weights while j<-i has no
intersections, only w(i->j) is recorded. This is determinded in code `if(u>0)`.
In this way, the undirected graph is symmetrized by halfing the weight
in code `weights(r, 2) = u/(2.0*ncol - u)/2`.
}
\details{
Author: Chen Hao, Date: 25/09/2015
Copied from https://github.com/JinmiaoChenLab/Rphenograph/blob/master/src/jaccard_coeff.cpp
for use in briDiscovr on 200611 under the Artistic-2.0 license
}
