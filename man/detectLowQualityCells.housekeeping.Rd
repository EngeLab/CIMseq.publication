% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{detectLowQualityCells.housekeeping}
\alias{detectLowQualityCells.housekeeping}
\title{detectLowQualityCells.housekeeping}
\arguments{
\item{counts}{data.frame; A data frame with counts data with gene names as
rownames and sample names as colnames.}

\item{geneName}{character; The gene name to use for the quantile cutoff. This
must be present in the rownames of the counts argument. Default is ACTB.}

\item{quantileCut}{numeric; This indicates probability at which the quantile
cutoff will be calculated using the normal distribution. Default = 0.01.}
}
\value{
A logical vector with length = ncol(counts) that is TRUE when the
counts data.frame column contains a sample with meeting the criteria specified
by the arguments.
}
\description{
It is often the case that some samples from sequencing experiments are of
low quality, in many cases due to issues during the sample preperation stage.
Due to the fact that these samples represent a high level of technical noise,
it is often desirable to remove these before downstream analysis which is
facilitated by this function. The function achieves this by utilizing a house
keeping gene and assuming its log2 expression to be normally distributed.
We then detect samples where the probability of the expression for the house
keeping gene in that sample is greater than the quantile.cut argument.
}
\examples{
set.seed(8292)
x <- runif(2e4)
y <- runif(2e4, 1.5, 100)
names <- paste0(letters, 1:2e4)
counts <- data.frame(a = x, b = y, c = y, row.names = names)
detectLowQualityCells.housekeeping(counts, geneName = "a1")

}
\author{
Jason Serviss
}
