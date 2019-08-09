# #context("permutationFunctions")
# 
# library(tidyverse)
# 
# ##run test permuteGroupMeans
# test_that("check that permuteGroupMeans outputs the expected result", {
#   
#   ###TEST1####
#   #prepare normal input data
#   means <- matrix(rep(11:13, each = 5), ncol = 3, nrow = 5)
#   
#   #setup expected data
#   nc <- ncol(means)
#   nr <- nrow(means)
#   
#   #run function
#   output <- permuteGroupMeans(means)
#   
#   #test
#   expect_equivalent(ncol(output), nc)
#   expect_equivalent(nrow(output), nr)
#   expect_equivalent(class(output), "matrix")
#   expect_false(identical(output, means))
# })
# 
# ##run test permuteSwarm
# test_that("check that permuteSwarm outputs the expected result", {
#   
#   ###TEST1####
#   #prepare normal input data
#   counts <- sp.scRNAseq::testCounts
#   sng <- str_detect(colnames(counts), "^s")
#   spCountsMul <- spCounts(counts[, !sng], matrix(NA, ncol = length(which(!sng))))
#   spUnsupervised <- sp.scRNAseq::testUns
#   iter <- 2
#   
#   #setup expected data
#   length <- l * nPerms * nrow(means)
#   value <-
#   
#   #run function
#   output <- permuteSwarm(
#     spCountsMul = spCountsMul,
#     spUnsupervised = spUnsupervised,
#     distFun = sp.scRNAseq:::distToSlice,
#     maxiter = 2,
#     swarmsize = 5,
#     cores = 1,
#     norm = TRUE,
#     iter = iter
#   )
#   
#   #test
#   expect_equivalent(length(output), length(1:iter))
# })
