% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{moveGenesToRownames}
\alias{moveGenesToRownames}
\title{moveGenesToRownames}
\arguments{
\item{counts}{data.frame; A data frame with counts data.}
}
\value{
The counts data.frame is returned with the first column moved to
the rownames of the data.frame.
}
\description{
Moves the first column of the counts data.frame to rownames and removes
the old column.
}
\examples{

counts <- data.frame(LETTERS, a = runif(26, 1, 100))
moveGenesToRownames(counts)

}
\author{
Jason Serviss
}
