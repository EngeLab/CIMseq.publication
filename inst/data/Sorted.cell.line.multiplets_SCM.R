#run from package root
#source('./inst/rawData/countsSorted2/Sorted.multiplets.R')

packages <- c("dplyr", "stringr", "CIMseq.publication")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

SCM <- function(save = TRUE) {
  projectName <- "Sorted.cell.line.multiplets_SCM"
  shortName <- "SCM"
  cat(paste0('Processing ', projectName, '\n'))

  m.link <- "https://drive.google.com/open?id=1XhVB0c9wi8BVKiDqhgp3sm4YNPnjD5im"
  Meta <- getMetadataViaLink(m.link)
  if("Missing" %in% colnames(Meta)) {
    Meta <- Meta %>%
      filter(is.na(Missing) | Missing == FALSE) %>%
      select(-Missing)
  }
  c.link <- 'https://drive.google.com/open?id=1jrfycO_evRcBPsk0cs7wYgF5L_7QG-eI'
  countData <- getCountsDataViaLinks(c.link)

  #move genes to rownames
  countData <- moveGenesToRownames(countData)

  #remove .htseq suffix
  if(any(grepl("htseq", colnames(countData)))) countData <- removeHTSEQsuffix(countData)

  #label singlets and multiplets ids should include SINGLET samples only
  singlets <- Meta[Meta$cellNumber == "Singlet", "sample"][[1]]
  singlets <- gsub("^..(.*)", "\\1", singlets)
  counts <- labelSingletsAndMultiplets(countData, singlets)

  #extract ERCC
  ercc <- detectERCCreads(counts)
  countsERCC <- counts[ercc, ]
  counts <- counts[!ercc, ]

  #remove non-genes
  counts <- counts[!detectNonGenes(counts), ]
  namesPreFilter <- colnames(counts)
  
  #filter data
  sng <- grepl("^s", colnames(counts))
  data.s <- filterCountsData(
    counts[, sng], countsERCC[, sng], geneMinCount = 0, cellMinCount = 4e4, 
    geneName = "ACTB", quantileCut.hk = 0.01, quantileCut.ercc = 0.95
  )
  data.m <- filterCountsData(
    counts[, !sng], countsERCC[, !sng], geneMinCount = 0, cellMinCount = 4.5e4, 
    geneName = "ACTB", quantileCut.hk = 0.01, quantileCut.ercc = 0.95
  )
  genes <- intersect(rownames(data.s[[1]]), rownames(data.m[[1]]))
  samples <- c(colnames(data.s[[1]]), colnames(data.m[[1]]))
  #add filtered column to Meta
  Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
    sample %in% samples, FALSE, TRUE
  ))

  #check all count samples in meta and vice versa
  c1 <- all(!Meta$sample %in% namesPreFilter)
  c2 <- all(!samples %in% Meta$sample)
  if(c1 & c2) {
    stop("all counts data not present in meta data")
  }

  #rename
  Counts <- cbind(data.s[[1]][genes, ], data.m[[1]][genes, ])
  CountsERCC <- cbind(data.s[[2]], data.m[[2]])

  #save .rda
  if(save) saveRDA(projectName, Counts, CountsERCC, Meta)
}
