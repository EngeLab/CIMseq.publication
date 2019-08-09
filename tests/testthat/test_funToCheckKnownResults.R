
context("funToCheckKnownResults")

##run test .tp
test_that("check that .tp outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- tibble::tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 5))),
    from = c("A", "B", rep("A", 3), rep("B", 2)),
    to = c("B", "C", "B", "C", "D", "C", "D")
  ) %>%
  tidyr::unite(connections, from, to, sep = "-") %>%
  tidyr::nest(-multiplet)
  
  known <- tibble::tibble(
    multiplet = sort(c(LETTERS[1], rep(LETTERS[2], 6))),
    from = c("A", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "B", "C", "D", "C", "D", "D")
  ) %>%
  tidyr::unite(connections, from, to, sep = "-") %>%
  tidyr::nest(-multiplet)
  
  data <- dplyr::full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(1, 5)
  
  #run function
  output <- .tp(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test .fp
test_that("check that .fp outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- tibble::tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 5))),
    from = c("A", "B", rep("A", 3), rep("B", 2)),
    to = c("B", "C", "B", "C", "D", "C", "D")
  ) %>%
    tidyr::unite(connections, from, to, sep = "-") %>%
    tidyr::nest(-multiplet)
  
  known <- tibble::tibble(
    multiplet = sort(c(LETTERS[1], rep(LETTERS[2], 6))),
    from = c("A", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "B", "C", "D", "C", "D", "D")
  ) %>%
    tidyr::unite(connections, from, to, sep = "-") %>%
    tidyr::nest(-multiplet)
  
  data <- dplyr::full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(1, 0)
  
  #run function
  output <- .fp(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test .fn
test_that("check that .fn outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- tibble::tibble(
    multiplet = sort(c(rep(LETTERS[1], 2), rep(LETTERS[2], 5))),
    from = c("A", "B", rep("A", 3), rep("B", 2)),
    to = c("B", "C", "B", "C", "D", "C", "D")
  ) %>%
    tidyr::unite(connections, from, to, sep = "-") %>%
    tidyr::nest(-multiplet)
  
  known <- tibble::tibble(
    multiplet = sort(c(LETTERS[1], rep(LETTERS[2], 6))),
    from = c("A", rep("A", 3), rep("B", 2), "C"),
    to = c("B", "B", "C", "D", "C", "D", "D")
  ) %>%
    tidyr::unite(connections, from, to, sep = "-") %>%
    tidyr::nest(-multiplet)
  
  data <- dplyr::full_join(detected, known, by = "multiplet")
  
  #setup expected data
  expected <- c(0, 1)
  
  #run function
  output <- .fn(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test .tn
test_that("check that .tn outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  detected <- CIMseq::getEdgesForMultiplet(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 
    rownames(getData(CIMseqSwarm_test, "fractions"))
  ) %>%
    tidyr::unite(connections, from, to, sep = "-") %>%
    dplyr::distinct() %>%
    tidyr::nest(-sample)
  
  known <- CIMseq::getEdgesForMultiplet(
    CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test, 
    rownames(getData(CIMseqSwarm_test, "fractions"))
  ) %>%
    tidyr::unite(connections, from, to, sep = "-") %>%
    dplyr::distinct() %>%
    tidyr::nest(-sample)
  
  data <- dplyr::full_join(detected, known, by = "sample")
  
  #setup expected data
  expected <- rep(2, 3)
  
  #run function
  output <- .tn(data, CIMseqSwarm_test, CIMseqSwarm_test)
  
  #test
  expect_equivalent(expected, output)
})

##run test .getCellTypes
test_that("check that .getCellTypes outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  data <- tibble::tibble(cellTypes = c("A1-B1", "A1-C1-D1", "A1-S1-F1-V1"))
  
  #setup expected data
  expected <- c("A1", "B1", "C1", "D1", "F1", "S1", "V1")
  
  #run function
  output <- .getCellTypes(data)
  
  #test
  expect_equivalent(expected, output)
})

##run test setupPlate
test_that("check that setupPlate outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  data <- tibble::tibble(
    row = c("A", "B"),
    column = 1:2,
    cellNumber = c("Multiplet", "Singlet"),
    sample = rep("test1", 2),
    cellTypes = rep("A1-B1", 2)
  )
  
  #setup expected data
  expected <- as.matrix(data.frame(A1 = 0.5, B1 = 0.5, row.names = "test1"))
  
  #run function
  output <- getData(setupPlate(data), "fractions")
  
  #test
  expect_equivalent(expected, output)
})

##run test viewAsPlate
test_that("check that viewAsPlate outputs the expected result", {
  
  ###TEST1####
  #prepare normal input data
  data <- structure(list(row = c("A", "B", "C", "D", "E"), column = c("01",
  "01", "01", "01", "01"), multipletComposition = c("HCT116-HOS", 
  "HCT116-HOS", "HCT116-HOS", "HCT116-HOS", "HCT116-HOS")), row.names = c(NA,
  -5L), class = c("tbl_df", "tbl", "data.frame"), .Names = c("row",
  "column", "multipletComposition"))
  
  #setup expected data
  expected <- tibble::tibble(
    row = LETTERS[1:5],
    `01` = rep("HCT116-HOS", 5)
  )
  
  #run function
  output <- viewAsPlate(data)
  
  #test
  expect_equivalent(expected, output)
})
