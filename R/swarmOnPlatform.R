#' runSwarmMultiprocess
#'
#'
#' @name runSwarmMultiprocess
#' @rdname runSwarmMultiprocess
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param swarmInit See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param maxiter See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param nSyntheticMultiplets See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param eps.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param maxit.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param currPath The analysis directory.
#' @param ... additional arguments to pass on
#' @return Saves swarm output to data folder in current directory.
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @import future.apply
#' @import future

runSwarmMultiprocess <- function(
  singlets, multiplets, 
  swarmInit = NULL, 
  maxiter = 100,
  nSyntheticMultiplets = 400,
  eps.stagnate = 1,
  maxit.stagnate = 5,
  currPath = getwd(),
  ...
){
  if(is.null(swarmInit)) swarmInit <- swarmInit(singlets, 2)
  plan(sequential)
  print(paste0("Starting deconvolution at ", Sys.time()))
  sObj <- CIMseqSwarm(
    singlets, multiplets, maxiter = maxiter, swarmsize = ncol(swarmInit), 
    nSyntheticMultiplets = nSyntheticMultiplets, swarmInit = swarmInit,
    psoControl = list(eps.stagnate = eps.stagnate, maxit.stagnate = maxit.stagnate)
  )
  print(paste0("Finished deconvolution at ", Sys.time()))
  if(!dir.exists('data')) dir.create('data')
  save(sObj, file = file.path(currPath, "data/sObj.rda"))
}

# .runSwarmBatch <- function(
#   singlets, multiplets, 
#   swarmInit = NULL, 
#   maxiter = 100,
#   nSyntheticMultiplets = 400,
#   eps.stagnate = 1,
#   maxit.stagnate = 5,
#   currPath = getwd(),
#   time = "24:00:00",
#   ...
# ){
#   
#   #requires future.batchtools package
#   if(is.null(swarmInit)) swarmInit <- swarmInit(singlets, 2)
#   options(future.wait.interval = 10000.0)
#   options(future.wait.timeout = 1e9)
#   plan(
#     future.batchtools::batchtools_slurm,
#     template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
#     resources = list(
#       account = "snic2018-8-151", partition = "core", ntasks = 1L,
#       time = time, jobname = "CIMseq",
#       modules = "R_packages/3.5.0", R = "R/3.5.0", 
#       log.file = file.path(currPath, "logs/slurm.txt")
#     ),
#     workers = 100
#   )
#   
#   #run deconvolution
#   print(paste0("Starting deconvolution at ", Sys.time()))
#   sObj <- CIMseqSwarm(
#     singlets, multiplets, maxiter = 100, swarmsize = 500, 
#     nSyntheticMultiplets = 400
#   )
#   print(paste0("Finished deconvolution at ", Sys.time()))
#   if(!dir.exists('data')) dir.create('data')
#   save(sObj, file = file.path(currPath, "data/sObj.rda"))
# }

#' runSwarmUppmax
#'
#'
#' @name runSwarmUppmax
#' @rdname runSwarmUppmax
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param swarmInit See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param maxiter See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param nSyntheticMultiplets See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param eps.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param maxit.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param currPath The analysis directory.
#' @param args List of arguments passed to the script. Should have sample in
#'  \code{args[1]}.
#' @param ... additional arguments to pass on
#' @return Saves swarm output to data folder in current directory.
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @import future.apply
#' @import future

runSwarmUppmax <- function(
  singlets, multiplets, 
  swarmInit = NULL, 
  maxiter = 100,
  nSyntheticMultiplets = 400,
  eps.stagnate = 1,
  maxit.stagnate = 5,
  e = 0.001,
  currPath = getwd(),
  args,
  ...
){
  plan(sequential)
  sample <- as.character(args[1])
  print(paste0("Running ", sample))
  multiplets.2 <- CIMseqMultiplets(
    getData(multiplets, "counts")[, sample, drop = FALSE],
    getData(multiplets, "counts.ercc")[, sample, drop = FALSE],
    getData(multiplets, "features")
  )
  
  #run deconvolution
  if(is.null(swarmInit)) swarmInit <- swarmInit(singlets, 2)
  print(paste0("Starting deconvolution at ", Sys.time()))
  sObj <- CIMseqSwarm(
    singlets, multiplets.2, maxiter = maxiter, swarmsize = ncol(swarmInit), 
    nSyntheticMultiplets = nSyntheticMultiplets, swarmInit = swarmInit, e = e,
    psoControl = list(eps.stagnate = eps.stagnate, maxit.stagnate = maxit.stagnate)
  )
  print(paste0("Finished deconvolution at ", Sys.time()))
  save(sObj, file = file.path(currPath, paste0("tmp/sObj_", sample, "_uppmax.rda")))
}