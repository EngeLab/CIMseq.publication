#' 
#' #' permuteGroupMeans
#' #'
#' #' Permutes the matrix in the groupMeans slot of the spUnsupervised object that
#' #' is input. To be called repeatedly to generate multiple permutations and used
#' #' for evaluating the performance of the spSwarm method. Permutation is per row
#' #' and, thus, the values in each row are randomly shuffled.
#' #'
#' #' @name permuteGroupMeans
#' #' @rdname permuteGroupMeans
#' #' @aliases permuteGroupMeans
#' #' @param spUnsupervised An spUnsupervised object.
#' #' @param ... additional arguments to pass on
#' #' @return Outputs one matrix of permuted group means with dims equal to the
#' #' origional group means matrix.
#' #' @author Jason T. Serviss
#' #'
#' #'
#' #'
#' NULL
#' #' @export
#' 
#' permuteGroupMeans <- function(
#'   groupMeans
#' ){
#'   t(apply(groupMeans, 1, sample))
#' }
#' 
#' #' permuteSwarm
#' #'
#' #' Write something
#' #'
#' #' Description
#' #'
#' #' @name permuteSwarm
#' #' @rdname permuteSwarm
#' #' @aliases permuteSwarm
#' #' @param spCountsMul An spCounts object containing multiplets.
#' #' @param spUnsupervised An spUnsupervised object.
#' #' @param distFun The distance function used to calculate the cost. Either the
#' #'    name of a custom function in the local environment or one of the included
#' #'    functions, i.e. \code{distToSlice, distToSliceNorm, distToSliceTop,
#' #'    distToSliceEuclid, distToSlicePearson, bic}.
#' #' @param maxiter pySwarm argument indicating maximum optimization iterations.
#' #' @param swarmsize pySwarm argument indicating the number of swarm particals.
#' #' @param cores The number of cores to be used while running spSwarm.
#' #' @param norm Logical indicating if the sum of fractions should equal 1.
#' #' @param iter The number of permutations to perform.
#' #' @param cellNumbers Tibble; Output from \code{estimateCells} function.
#' #' @param e Numeric; The epsilon value for the .complexityPenilty unit.
#' #' @param ... additional arguments to pass on
#' #' @return ?.
#' #' @author Jason T. Serviss
#' #' @keywords permuteSwarm
#' #' @examples
#' #'
#' #' #use demo data
#' #'
#' #'
#' #' #run function
#' #'
#' #'
#' NULL
#' 
#' #' @rdname permuteSwarm
#' #' @export
#' 
#' setGeneric("permuteSwarm", function(
#'     spCountsMul,
#'     ...
#' ){
#'     standardGeneric("permuteSwarm")
#' })
#' 
#' #' @rdname permuteSwarm
#' #' @export
#' #' @import CIMseq
#' #' @importFrom dplyr pull
#' #' @importFrom parallel mclapply
#' 
#' setMethod("permuteSwarm", "CIMseqMultiplets", function(
#'   spCountsMul,
#'   spUnsupervised,
#'   distFun = sp.scRNAseq:::dtsnCellNum,
#'   maxiter = 10,
#'   swarmsize = 50,
#'   cores = 1,
#'   norm = TRUE,
#'   iter,
#'   cellNumbers,
#'   e = 0.0025,
#'   ...
#' ){
#'   distFun <- match.fun(distFun)
#'   geneIdx <- getData(spUnsupervised, "selectInd")
#'   groupMeans <- getData(spUnsupervised, "groupMeans")[geneIdx, ]
#'   multiplets <- getData(spCountsMul, "counts.cpm")[geneIdx, ]
#'   
#'   fractions <- rep(1.0 / ncol(groupMeans), ncol(groupMeans))
#'   to <- if(ncol(multiplets) == 1) {to <- 1} else {to <- ncol(multiplets)}
#'   control <- list(maxit = maxiter, s = swarmsize)
#'   matchIdx <- match(colnames(multiplets), pull(cellNumbers, sampleName))
#'   cellNumbers <- pull(cellNumbers, cellNumberMedian)[matchIdx]
#'   
#'   permData <- .runPermutations(
#'     iter,
#'     multiplets,
#'     groupMeans,
#'     distFun = distFun,
#'     fractions,
#'     to,
#'     control,
#'     cores = cores,
#'     cellNumbers,
#'     e,
#'     ...
#'   )
#'   
#'   return(permData)
#' })
#' 
#' .runPermutations <- function(
#'   iter,
#'   multiplets,
#'   groupMeans,
#'   distFun = distFun,
#'   fractions,
#'   to,
#'   control,
#'   cores,
#'   cellNumbers,
#'   e,
#'   ...
#' ){
#'   permData <- list()
#'   
#'   for(y in 1:iter) {
#'     set.seed(y)
#'     permData[[y]] <- mclapply(1:to, function(i, ...) {
#'       sp.scRNAseq:::.optim.fn(
#'         i,
#'         fractions,
#'         distFun,
#'         permuteGroupMeans(groupMeans),
#'         control,
#'         multiplets,
#'         cellNumbers,
#'         e,
#'         ...
#'       )
#'     }, mc.cores = cores)
#'   }
#'   return(permData)
#' }
#' 
#' #' spSwarmPermutationData
#' #'
#' #' A function to put the results from the \code{permuteSwarm} function into
#' #' spSwarm objects.
#' #'
#' #' @name spSwarmPermutationData
#' #' @rdname spSwarmPermutationData
#' #' @aliases spSwarmPermutationData
#' #' @param sObj An spSwarm object where the "real" data is included.
#' #' @param permData Output from the permuteSwarm function.
#' #' @param ... additional arguments to pass on
#' #' @return A list of spSwarm objects.
#' #' @author Jason T. Serviss
#' #'
#' #'
#' #'
#' NULL
#' #' @export
#' #' @import CIMseq
#' #' @importFrom purrr map
#' 
#' spSwarmPermutationData <- function(sObj, permData) {
#'   map(permData, .onespSwarmObj, sObj)
#' }
#' 
#' .onespSwarmObj <- function(data, sObj) {
#'   output <- data.frame(t(sapply(data, function(j) j[[1]])))
#'   cost <- sapply(data, function(j) j[[2]])
#'   counts <- t(sapply(data, function(j) j[[3]]))
#'   convergence <- sapply(data, function(j) j[[4]])
#'   convergenceKey <- c(
#'     "Maximal number of function evaluations reached." = 1,
#'     "Maximal number of iterations reached." = 2,
#'     "Maximal number of restarts reached." = 3,
#'     "Maximal number of iterations without improvement reached." = 4
#'   )
#'   convergence <- names(convergenceKey)[match(convergence, convergenceKey)]
#'   
#'   
#'   #normalize swarm output
#'   output <- output * 1/rowSums(output)
#'   
#'   #add names
#'   spSwarm <- getData(sObj, "spSwarm")
#'   colnames(output) <- colnames(spSwarm)
#'   rownames(output) <- rownames(spSwarm)
#'     
#'   #create object
#'   new("spSwarm",
#'       spSwarm = output,
#'       costs = cost,
#'       convergence = convergence,
#'       stats = list(),
#'       arguments = getData(sObj, "arguments")
#'   )
#' }
#' 
#' #' tidyPermutationData
#' #'
#' #' A function to wrangle the results from the \code{permuteSwarm} function.
#' #'
#' #' @name tidyPermutationData
#' #' @rdname tidyPermutationData
#' #' @aliases tidyPermutationData
#' #' @param sObj An spSwarm object where the "real" data is included.
#' #' @param permData Output from the permuteSwarm function.
#' #' @param ... additional arguments to pass on
#' #' @return A tibble with the real and permuted costs and fractions.
#' #' @author Jason T. Serviss
#' #'
#' #'
#' #'
#' NULL
#' #' @export
#' #' @import CIMseq
#' #' @importFrom tibble tibble as_tibble
#' #' @importFrom dplyr mutate bind_rows rename select group_by full_join
#' #' @importFrom purrr map map_dbl simplify_all
#' #' @importFrom tidyr nest
#' 
#' tidyPermutationData <- function(sObj, permData) {
#'   #setup real data
#'   spSwarm <- getData(sObj, "spSwarm")
#'   real <- tibble(
#'     multiplet = rownames(spSwarm),
#'     cost = getData(sObj, "costs")
#'   ) %>%
#'   mutate(fracs = map(seq_len(nrow(spSwarm)), function(i) spSwarm[i, ]))
#'   
#'   #get cost and fractions (normalized) and ? from permData
#'   #cost (called "value" in permData object)
#'   costs <- map(permData, function(x) {
#'     map_dbl(x, function(y) {
#'       y[['value']]
#'     })
#'   }) %>%
#'   map(., as_tibble) %>%
#'   bind_rows() %>%
#'   rename(cost = value) %>%
#'   mutate(multiplet = rep(rownames(spSwarm), length(permData))) %>%
#'   group_by(multiplet) %>%
#'   nest(cost, .key = "permCosts")
#'   
#'   #fractions
#'   fracs <- map(permData, function(x) {
#'     map(x, function(y) {
#'       y[['par']]
#'     })
#'   }) %>%
#'   unlist() %>%
#'   matrix(., ncol = 3) %>%
#'   `*` (1/rowSums(.)) %>%
#'   as_tibble() %>%
#'   mutate(multiplet = rep(rownames(spSwarm), length(permData))) %>%
#'   group_by(multiplet) %>%
#'   nest(-multiplet, .key = "permFracs")
#'   
#'   #combine all
#'   real %>%
#'     full_join(costs) %>%
#'     full_join(fracs)
#' }
#' 
#' #' permCostPlot
#' #'
#' #' A function to plot the permuted vs real costs.
#' #'
#' #' @name permCostPlot
#' #' @rdname permCostPlot
#' #' @aliases permCostPlot
#' #' @param data Output from the tidyPermutationData function.
#' #' @param ... additional arguments to pass on
#' #' @return A tibble with the real and permuted costs and fractions.
#' #' @author Jason T. Serviss
#' #'
#' #'
#' #'
#' NULL
#' #' @export
#' #' @import CIMseq
#' #' @import ggplot2
#' #' @importFrom ggthemes theme_few
#' #' @importFrom dplyr select
#' #' @importFrom tidyr unnest
#' 
#' permCostPlot <- function(data, ...) {
#'   #setup "real" data
#'   real <- select(data, multiplet, cost)
#'   
#'   #process permuted data and plot
#'   p <- data %>%
#'   select(multiplet, permCosts) %>%
#'   unnest() %>%
#'   ggplot(data = .,) +
#'   geom_boxplot(aes(x = multiplet, y = cost)) +
#'   geom_point(data = real, aes(x = multiplet, y = cost), colour = "red") +
#'   theme_few() +
#'   theme(axis.text.x = element_text(angle = 90))
#'   
#'   p
#'   return(p)
#' }
#' 
#' #' calculateCostP
#' #'
#' #' A function to plot the permuted vs real costs.
#' #'
#' #' @name calculateCostP
#' #' @rdname calculateCostP
#' #' @aliases calculateCostP
#' #' @param data Output from the tidyPermutationData function.
#' #' @param ... additional arguments to pass on
#' #' @return A tibble with the real and permuted costs and fractions.
#' #' @author Jason T. Serviss
#' #'
#' #'
#' #'
#' NULL
#' #' @export
#' #' @import CIMseq
#' #' @import ggplot2
#' #' @importFrom purrr pmap simplify
#' #' @importFrom dplyr if_else pull "%>%"
#' #' @importFrom tibble tibble
#' 
#' calculateCostP <- function(data) {
#'   data %>%
#'   {pmap(list(.$cost, .$permCosts), function(cost, permCost) {
#'     permCost <- as.numeric(permCost[[1]])
#'     l <- length(permCost)
#'     gt <- sum(permCost < cost) / l
#'     if_else(gt == 0, 10^-log10(l), gt)
#'   })} %>%
#'   simplify() %>%
#'   tibble(
#'     multiplet = pull(data, multiplet),
#'     pValue = .,
#'   )
#' }
