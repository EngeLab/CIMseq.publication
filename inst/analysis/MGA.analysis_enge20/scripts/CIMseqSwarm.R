#ARGS
args <- commandArgs(TRUE)

#PACKAGES
packages <- c("tidyverse", "future")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
if(file.exists(file.path(currPath, 'data/CIMseqData.rda'))) {
  load(file.path(currPath, 'data/CIMseqData.rda'))
}

if(!is.na(args[1])) {
  library(CIMseq, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  library(CIMseq.data, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  library(CIMseq.testing, lib.loc = "/home/jason/R/x86_64-redhat-linux-gnu-library/3.5")
  out <- if(is.na(args[2])) {out <- 'tmp'} else {out <- as.character(args[2])}
  
  #gives 1370 swarm members
  baseSeed <- 43892
  init <- cbind(
      swarmInit(cObjSng, 2, null.weight = 0.5, seed = baseSeed),
      swarmInit(cObjSng, 3, null.weight = 0.5, seed = baseSeed)
  )
  options(future.globals.maxSize = Inf)
  runSwarmUppmax(
    cObjSng, cObjMul, swarmInit = init, maxiter = 100, swarmsize = ncol(init),
    nSyntheticMultiplets = 2000, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = out, args = args
  )
  
} else {
  library(CIMseq)
  library(CIMseq.publication)
  
  #gives 1370 swarm members
  baseSeed <- 43892
  init <- swarmInit(cObjSng, 2, null.weight = 0.5, seed = baseSeed)
  
  options(future.globals.maxSize = Inf)
  runSwarmMultiprocess(
    cObjSng, cObjMul, swarmInit = init, maxiter = 2, swarmsize = ncol(init),
    nSyntheticMultiplets = 5, eps.stagnate = 1, maxit.stagnate = 5,
    currPath = getwd()
  )
}

writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqSwarm.txt"))
print("finished")