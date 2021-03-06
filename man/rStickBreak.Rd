% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infinite.R
\name{rStickBreak}
\alias{rStickBreak}
\title{Calculate random proportions from a dirchlet process based on stick breaking procedure}
\usage{
rStickBreak(n, alpha = 1)
}
\arguments{
\item{n}{number of numbers to generate}

\item{alpha}{smaller alpha means on average more weight concentrated earlier}
}
\value{
a vector of n numbers adding to 1
}
\description{
Calculate random proportions from a dirchlet process based on stick breaking procedure
}
\examples{
rStickBreak(10)
rStickBreak(10,3)
}

