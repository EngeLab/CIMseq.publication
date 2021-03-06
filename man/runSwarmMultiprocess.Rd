% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swarmOnPlatform.R
\name{runSwarmMultiprocess}
\alias{runSwarmMultiprocess}
\title{runSwarmMultiprocess}
\arguments{
\item{singlets}{CIMseqSinglets; A CIMseqSinglets object.}

\item{multiplets}{CIMseqMultiplets; A CIMseqMultiplets object.}

\item{swarmInit}{See \code{\link[CIMseq]{CIMseqSwarm}}.}

\item{maxiter}{See \code{\link[CIMseq]{CIMseqSwarm}}.}

\item{nSyntheticMultiplets}{See \code{\link[CIMseq]{CIMseqSwarm}}.}

\item{eps.stagnate}{See \code{\link[CIMseq]{pso.2.0}}.}

\item{maxit.stagnate}{See \code{\link[CIMseq]{pso.2.0}}.}

\item{currPath}{The analysis directory.}

\item{...}{additional arguments to pass on}
}
\value{
Saves swarm output to data folder in current directory.
}
\description{
runSwarmMultiprocess
}
\author{
Jason T. Serviss
}
