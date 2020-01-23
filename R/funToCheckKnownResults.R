
#' checkResults
#'
#' Untested for self connections.
#'
#' @name checkResults
#' @rdname checkResults
#' @aliases checkResults
#' @param swarm CIMseqSwarm; A CIMseqSwarm object.
#' @param known A CIMseqSwarm object of the known result. The results should be
#' represented in the fractions slot according to the edge.cutoff provided.
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @importFrom tidyr unite
#' @importFrom dplyr group_by summarize "%>%" mutate rename everything select
#' @importFrom purrr map2_int pmap_int map_chr
#' @importFrom tibble add_column as_tibble
#' @importFrom utils combn

checkResults <- function(
  swarm, known, singlets, multiplets, ...
){
  connections <- from <- to <- tp <- fn <- tn <- fp <- data.x <- NULL
  data.y <- cells <- data.detected <- data.expected <- NULL
  
  detected <- CIMseq::getCellsForMultiplet(
    swarm, singlets, multiplets
  ) %>%
    group_by(sample) %>%
    summarize(data = list(cells))
  
  expected <- CIMseq::getCellsForMultiplet(
    known, singlets, multiplets
  ) %>%
    group_by(sample) %>%
    summarize(data = list(cells))
  
  #multiplets with 0 of 1 (self) connections will not be included in output from 
  # getCellsForMultiplet.
  
  full.data <- full_join(detected, expected, by = "sample")
  
  #check here that all muultiplets are included in detected and expected
  if(!all(colnames(getData(multiplets, "counts")) %in% pull(full.data, sample))) {
    stop("detected missing samples. check checkResults FUN")
  }
  
  full.data %>%
    mutate(
      tp = .tp(.),
      fp = .fp(.),
      fn = .fn(.),
      tn = .tn(., swarm, known),
      TPR = .TPR(tp, fn),
      TNR = .TNR(tn, fp),
      ACC = .ACC(tp, tn, fp, fn)
    ) %>%
    rename(data.detected = data.x, data.expected = data.y) %>%
    mutate(data.detected = map_chr(data.detected, ~paste(.x, collapse = ", "))) %>%
    mutate(data.expected = map_chr(data.expected, ~paste(.x, collapse = ", "))) %>%
    select(sample, data.expected, data.detected, everything())
}

#calculates the true positives
.tp <- function(data) {
  data %>%
  {map2_int(.$data.x, .$data.y, function(detected, expected) {
    sum(unlist(expected) %in% unlist(detected))
  })}
}

#calculates the false positives
.fp <- function(data) {
  data %>%
  {map2_int(.$data.x, .$data.y, function(detected, expected) {
    sum(!unlist(detected) %in% unlist(expected))
  })}
}

#calculates the false negatives
.fn <- function(data) {
  data %>%
  {map2_int(.$data.x, .$data.y, function(detected, expected) {
    sum(!unlist(expected) %in% unlist(detected))
  })}
}

#calculates the true negatives
#in order to know which combinations are possible and which are missing from
# both expected and detected, the .possibleCombs function is used.
.tn <- function(data, swarm, known) {
  possibleCombs <- .possibleCombs(swarm, known)

  data %>%
  add_column(possibleCombs = list(possibleCombs)) %>%
  {pmap_int(
    list(.$data.x, .$data.y, .$possibleCombs),
    function(detected, expected, possibleCombs) {
      bool1 <- !unlist(possibleCombs) %in% unlist(expected)
      bool2 <- !unlist(possibleCombs) %in% unlist(detected)
      sum(bool1 & bool2)
    }
  )}
}

#Calculates all possible connections with all cell types
#Helper for .tn
# .possibleCombs <- function(swarm, known) {
#   connections <- V1 <- V2 <- NULL
#   ctKnown <- colnames(getData(known, "fractions"))
#   ctDetected <- colnames(getData(swarm, "fractions"))
# 
#   c(ctKnown, ctDetected) %>%
#     unique() %>%
#     combn(., 2) %>%
#     t() %>%
#     as_tibble() %>%
#     unite(connections, V1, V2, sep = "-")
# }

.possibleCombs <- function(swarm, known) {
  connections <- V1 <- V2 <- NULL
  ctKnown <- colnames(getData(known, "fractions"))
  ctDetected <- colnames(getData(swarm, "fractions"))
  unique(c(ctKnown, ctDetected))
}

#Calculates the true positive rate (sensitivity)
.TPR <- function(tp, fn) {
  tp / (tp + fn)
}

#Calculates the true negative rate (specificity)
.TNR <- function(tn, fp) {
  tn / (tn + fp)
}

#Calculates the accuracy
.ACC <- function(tp, tn, fp, fn) {
  (tp + tn) / (tp + fn + fp + tn)
}

#' printResults
#'
#' Takes the output from checkResults and prints it in a more verbose form.
#' (note: rewrite as general to take an arg for nested columns to expand)
#'
#' @name printResults
#' @rdname printResults
#' @aliases printResults
#' @param data Output from checkResults function.
#' @param addFractions Logical indicating if fractions should be added to the
#'  output.
#' @param spSwarm An spSwarm object. Required if \code{addFractions} is TRUE.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @importFrom dplyr pull "%>%" arrange bind_rows group_by summarize rename
#' @importFrom stringr str_split
#' @importFrom purrr map_dfr
#' @importFrom tibble add_column column_to_rownames rownames_to_column
#' @importFrom stats setNames
#' @importFrom tidyr unnest

printResults <- function(
  data, addFractions = TRUE, spSwarm = NULL, ...
){
  multiplet <- connections.x <- connections.y <- costFunction <- NULL
  cellsInWell <- tp <- ACC <- NULL
  #expand connections
  data.detected <- select(data, multiplet, data.detected) %>%
  .expandNested()

  data.expected <- select(data, multiplet, data.expected) %>%
  .expandNested()

  expanded <- select(data, -data.detected, -data.expected) %>%
  full_join(data.detected, by = "multiplet") %>%
  full_join(data.expected, by = "multiplet") %>%
  rename(data.detected = connections.x, data.expected = connections.y) %>%
  select(costFunction:cellsInWell, data.expected, data.detected, tp:ACC) %>%
  arrange(costFunction, multiplet)

  #add fractions
  if(addFractions) {
    return(.addFractions(expanded, spSwarm))
  } else {
    return(expanded)
  }
}

.expandNested <- function(contracted) {
  multiplet <- connections <- NULL
  contracted %>%
  unnest() %>%
  bind_rows() %>%
  group_by(multiplet) %>% #change in generalized solution
  summarize(connections = paste(connections, collapse = ", "))
}

.addFractions <- function(expanded, spSwarm) {
  fractions <- multiplet <- cellsInWell <- data.expected <- ACC <- NULL
  fracs <- getData(spSwarm, "spSwarm")
  names <- paste(colnames(fracs), collapse = ", ")
  formated <- paste("frac (", names, ")", sep = "")

  fracs %>%
  round(digits = 2) %>%
  rownames_to_column("multiplet") %>%
  as_tibble() %>%
  full_join(expanded, by = "multiplet") %>%
  unite(fractions, 2:(ncol(fracs) + 1), sep = ", ") %>%
  setNames(c(colnames(.)[1], formated, colnames(.)[3:ncol(.)])) %>%
  select(multiplet, cellsInWell, 2, data.expected:ACC) %>%
  arrange(multiplet)
}


#' resultsInPlate
#'
#' Takes the output from checkResults and setupPlate and displays results in
#' plate format.
#'
#' @name resultsInPlate
#' @rdname resultsInPlate
#' @aliases resultsInPlate
#' @param results Output from checkResults function.
#' @param plate Output from setup plate or a tibble including the columns \code{
#'  row, column, and multipletName} indicating the position of each multiplet in
#'  the plate.
#' @param var The variable to visualize. May be "tp, fp, tn, fn, TPR, TNR, ACC".
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss

NULL

#' @export
#' @import CIMseq
#' @importFrom dplyr select full_join "%>%" rename
#' @import ggplot2

resultsInPlate <- function(
  results,
  plate,
  var,
  ...
){
  multiplet <- column <- NULL
  results %>%
  select(multiplet, var) %>%
  full_join(plate, by = c("multiplet" = "multipletName")) %>%
  ggplot(
    data = .,
    aes(x = column, y = ordered(row, levels = rev(LETTERS[1:8])))
  ) +
  geom_tile(aes_string(fill = var)) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.position = "top"
  ) +
  guides(fill = guide_colourbar(barwidth = 13))
}

#' setupPlate
#'
#' Sets up the CIMseqSwarm object for the \code{known} argument to the
#' \code{checkResults} function.
#'
#' @name setupPlate
#' @rdname setupPlate
#' @aliases setupPlate
#' @param plateData tibble; Includes the columns \code{row, column,
#'  multipletName, multipletComposition, connections} where \code{row} and
#'  \code{column} includes a letter or numerical value, respectivley, indicating
#'  the index in the plate of each multiplet, \code{multipletName} corresponds
#'  to the name of the multiplet and matches the names in the subsequent
#'  detected results, \code{multipletComposition} is the cell types included in
#'  the multiplet seperated by a "-", and \code{connections} is a list for each
#'  multiplet containing a matrix with all possible combinations of 2 of the
#'  cell types included in the multiplet.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss
#'
#'
#'
NULL

#' @export
#' @import CIMseq
#' @importFrom dplyr pull "%>%" filter mutate select
#' @importFrom stringr str_split
#' @importFrom purrr map_dfr map map2
#' @importFrom tibble tibble add_column column_to_rownames
#' @importFrom utils combn
#' @importFrom methods new

setupPlate <- function(
  plateData,
  fill = NULL,
  ...
){
  cellNumber <- cellTypes <- connections <- fractions <- NULL
  
  #make CIMseqSwarm slot for CIMseqSwarm object
  u.classes <- pull(plateData, cellTypes) %>%
    str_split(., "-") %>%
    unlist() %>%
    unique() %>%
    sort()
  
  spSwarm <- plateData %>%
    select(sample, cellNumber, cellTypes) %>%
    filter(cellNumber == "Multiplet") %>%
    mutate(connections = str_split(cellTypes, "-")) %>%
    mutate(fractions = map2(connections, sample, function(c, s) {
      l <- length(unique(c))
      frac <- rep(0, length(u.classes))
      names(frac) <- u.classes
      if(is.null(fill)) frac[c] <- 1/l else frac[c] <- 1
      m <- matrix(frac, ncol = length(frac), dimnames= list(NULL, names(frac)))
      data.frame(sample = s, m, stringsAsFactors = FALSE)
    })) %>%
    select(fractions) %>%
    unnest(cols = c(fractions)) %>%
    as.data.frame() %>%
    column_to_rownames("sample") %>%
    as.matrix()

  #create spSwarm object
  new("CIMseqSwarm",
    fractions = as.matrix(spSwarm),
    costs = vector(mode = "numeric"),
    convergence = vector(mode = "character"),
    stats = tibble(),
    singletIdx = list(),
    arguments = tibble()
  )
}

.getCellTypes <- function(plateData) {
  cellTypes <- NULL
  plateData %>%
    pull(cellTypes) %>%
    str_split("-") %>%
    unlist() %>%
    unique() %>%
    sort()
}

#' viewAsPlate
#'
#' Takes the output from setupPlate and reformats it to a 96-well plate format.
#'
#' @name viewAsPlate
#' @rdname viewAsPlate
#' @aliases viewAsPlate
#' @param plate Output from setup plate.
#' @param ... additional arguments to pass on
#' @return data.frame
#' @author Jason T. Serviss

NULL

#' @export
#' @importFrom dplyr select "%>%"
#' @importFrom tidyr spread

viewAsPlate <- function(plate) {
  column <- multipletComposition <- NULL
  plate %>%
  select(row, column, multipletComposition) %>%
  spread(column, multipletComposition)
}
