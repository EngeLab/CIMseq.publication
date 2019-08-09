#' Reproduction of the CIMseq publication figures.
#'
#'
#' \tabular{ll}{ Package: \tab CIMseq.publication\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2016-02-28\cr License: \tab GPL-3\cr }
#'
#' @name CIMseq.publication-package
#' @aliases CIMseq.publication-package CIMseq.publication
#' @docType package
#' @author Author: Jason T. Serviss
#' @keywords package
#'
#' @import methods
#' @import ggplot2
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
