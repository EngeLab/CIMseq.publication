
#' syntheticSinglets
#'
#' This unit uses the negative binomial distribution to synthesize the
#' singlets.
#'
#'
#' @name syntheticSinglets
#' @aliases syntheticSinglets
#' @param nGenes Number of genes in generated synthetic data.
#' @param nCells Number of cells per cell type in generated synthetic data.
#' @param nCellTypes Number of cell types in generated synthetic data.
#' @param ... additional arguments to pass on
#' @return A matrix with synthetic counts.
#' @author Jason T. Serviss
#' @keywords internal
#' @examples
#'
#' synth <- syntheticSinglets(10, 10, 10)
#'
NULL
#' @export
#' @importFrom stats rnbinom runif

syntheticSinglets <- function(
  nGenes, nCells, nCellTypes, seed = 8732, ...
){
  synth <- sapply(1:nCellTypes, function(x) {
    set.seed(seed + x)
    rnbinom(nGenes * nCells, mu = 2^runif(nGenes, 0, 5), size = x)
  })
    
  singlets <- matrix(as.numeric(synth), nrow = nGenes)
    
  colnames(singlets) <- paste(
    sort(rep(LETTERS, nCells))[1:(nCellTypes * nCells)],
    1:nCells, sep = ""
  )
  as.matrix(singlets)
}

#' syntheticMultipletsFromCounts
#'
#' Generates one multiplet per combo using the singlet counts data input to the
#' function.
#'
#' @name syntheticMultipletsFromCounts
#' @rdname syntheticMultipletsFromCounts
#' @aliases syntheticMultipletsFromCounts
#' @param counts matrix; A matrix of singlet counts per million with cells as
#'  columns and genes as rows.
#' @param classes character; A character vector indicating the cell types of the
#'  cells in the counts matrix.
#' @param combos matrix; A matrix the combinations of cell types to be included
#'  in the multiplets. Typically generated with the \code{combn} function.
#' @param adjustment; numeric A named numeric vector with length = the number of
#'  unique cell types indicating the fraction of contribution for each cell
#'  type. Names should indicate the corresponding cell type.
#' @param seed numeric; Seed for random number initiation. When generating many
#'  multiplets by running the function multiple times in a loop, the seed should
#'  be set dynamically in the function calling the makeSyntheticData function to
#'  avoid all multiplets being identical.
#' @param ... Additional arguments to pass on.
#' @return A tibble with one gene per row and one synthesized multiplet per
#'  combo as columns.
#' @author Jason T. Serviss
#' @examples
#'
#' #use demo data
#'
#'
NULL
#' @export
#' @import CIMseq
#' @importFrom dplyr "%>%" full_join
#' @importFrom purrr map map_int reduce
#' @importFrom tibble as_tibble

syntheticMultipletsFromCounts <- function(
  counts, classes, fractions, seed = 87909023, ...
){
  set.seed(seed)
  out <- as.numeric(sampleSinglets(classes)) %>%
  subsetSinglets(counts, .) %>%
  adjustAccordingToFractions(fractions, .) %>%
  multipletSums()
  
  rownames(out) <- rownames(counts)
  out
}

#syntheticMultipletsFromCounts <- function(
#  counts,
#  classes,
#  combos,
#  adjustment,
#  seed = 87909023,
#  ...
#){
#  set.seed(seed)
#
#  #setup output structure
#  output <- data.frame(gene = rownames(counts), stringsAsFactors = FALSE)
#
#  combos %>%
#  as.data.frame(stringsAsFactors = FALSE) %>%
#  #for each combination in combos...
#  map(., function(x) {
#    #the x variable is a character vector including the names of the cell types
#    # to be included in the multiplet
#
#    #select one cell from the pool for each cell type
#    exCounts <- counts[, map_int(x, ~sample(which(classes == .x), 1))]
#   colnames(exCounts) <- x
#
#    #adjust the counts by multiplying with the corresponding values in the
#    # adjustment variable provided as an arg
#    adj <- adjustment[match(x, names(adjustment))]
#    adjusted <- round(t(t(exCounts) * adj))
#
#    #calculate the sum of counts for each gene and expand the gene names to the
#    # length of the corresponding counts. This provides the total pool of counts
#    # to sample from to generate the multiplet.
#
#    rs <- matrixStats::rowSums2(adjusted)
#
#    #Sample from the poisson distribution for each gene with lambda = rs
#    matrix(
#      rpois(n = length(rs), lambda = rs),
#      dimnames = list(rownames(exCounts), "multiplet")
#    ) %>%
#      as.data.frame() %>%
#      rownames_to_column("gene") %>%
#      setNames(c("gene", paste(sort(x), collapse = "-")))
#  }) %>%
#  #reformat as tibble for easy integration with output from other function calls
#  reduce(full_join, by = "gene") %>%
#  full_join(output, by = "gene") %>%
#  replace(is.na(.), 0) %>%
#  as_tibble()
#}



