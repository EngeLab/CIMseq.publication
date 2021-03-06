% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funToCheckKnownResults.R
\name{resultsInPlate}
\alias{resultsInPlate}
\title{resultsInPlate}
\arguments{
\item{results}{Output from checkResults function.}

\item{plate}{Output from setup plate or a tibble including the columns \code{
row, column, and multipletName} indicating the position of each multiplet in
the plate.}

\item{var}{The variable to visualize. May be "tp, fp, tn, fn, TPR, TNR, ACC".}

\item{...}{additional arguments to pass on}
}
\value{
data.frame
}
\description{
Takes the output from checkResults and setupPlate and displays results in
plate format.
}
\author{
Jason T. Serviss
}
