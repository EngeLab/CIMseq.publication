#run from package root
#source('inst/rawData/countsMgfp/Mouse.gut.architecture.R')

packages <- c("dplyr", "stringr", "CIMseq.publication")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

MGA <- function(save = TRUE) {
  projectName <- "Mouse.gut.architecture_MGA"
  shortName <- "MGA"
  cat(paste0('Processing ', projectName, '\n'))
  
  m.link <- 'https://drive.google.com/open?id=18eCfbpbt5lagGz3gnz3Rs_T-BXe-Gwsu'
  Meta <- getMetadataViaLink(m.link)
  if("Missing" %in% colnames(Meta)) {
    Meta <- Meta %>%
      filter(is.na(Missing) | Missing == FALSE) %>%
      select(-Missing)
  }
  
  prefix <- 'https://drive.google.com/file/d/'
  suffix <- '/view?usp=sharing'
  
  short.id <- c(
    '1y0CgZbSAK4TnVtozBPgVIpgTBeQZoEn3',
    '1FFedVrk1ajeV4U4AI2Da37gH2pzfnNWi',
    '1ztgojZRaMqb54tlLTiXyAMNy8Cso42On',
    '1t6rycfkZRHEX1ZGnQ_stwawWDZu9Y-HZ',
    '1P8nnel-jT8_6Qgrpg59Htmfo19-UDYXD',
    '1g8qi6MRMygspluBNDdcEIqUnkcHoh3oJ',
    '1d_EPDSIf5V7qsAvmFZOuWkjzxp0V4xiQ',
    '1f4T1z6B5-Yph2QK1o2bWkeFh4e0h3JjX',
    '12IdD7Z5WQ0qMxFBOUsO-CiftZj09kYOR',
    '1rNyi3KAuTrQf5kgaKVyWLSRd7fLyTxbH',
    '1k-wol_MoHp7kObRrhyVVJTvHClP6P84f',
    '1vsaSqzR_whbJffSgsVm5YuHKfiBHL_Ke',
    '1OVoKJGVZDMXRAX82vF-JPjBfIcQM3RpK',
    '1tktL-51WY30eHw6Fa17fpEo1SHjxGwWI'
  )
  
  c.link <- paste(prefix, short.id, suffix, sep = "")
  
  countData <- getCountsDataViaLinks(c.link)
  colnames(countData) <- renameMgfpSamples(colnames(countData))

  #move genes to rownames
  countData <- moveGenesToRownames(countData)

  #label singlets and multiplets ids should include SINGLET samples only
  singlets <- Meta[Meta$cellNumber == "Singlet", "sample"][[1]]
  singlets <- gsub("^..(.*)", "\\1", singlets)
  countData <- labelSingletsAndMultiplets(countData, singlets)

  #extract ERCC
  ercc <- detectERCCreads(countData)
  CountsERCC <- countData[ercc, ]
  Counts <- countData[!ercc, ]

  #remove non-genes
  Counts <- Counts[!detectNonGenes(Counts), ]

  #filter counts
  sng <- grepl("^s", colnames(Counts))
  data.s <- filterCountsData(
    Counts[, sng], CountsERCC[, sng], geneMinCount = 0, cellMinCount = 1.8e4, 
    geneName = "Actb", quantileCut.hk = 0.01, quantileCut.ercc = 0.99
  )
  data.m <- filterCountsData(
    Counts[, !sng], CountsERCC[, !sng], geneMinCount = 0, cellMinCount = 2.3e4, 
    geneName = "Actb", quantileCut.hk = 0.01, quantileCut.ercc = 0.99
  )
  genes <- intersect(rownames(data.s[[1]]), rownames(data.m[[1]]))
  samples <- c(colnames(data.s[[1]]), colnames(data.m[[1]]))
  
  #add filtered column to Meta
  Meta <- dplyr::mutate(Meta, filtered = dplyr::if_else(
    sample %in% samples, FALSE, TRUE
  ))

  #check all count samples in meta and vice versa
  c1 <- all(!Meta$sample %in% colnames(Counts))
  c2 <- all(!samples %in% Meta$sample)
  if(c1 & c2) {
    stop("all counts data not present in meta data")
  }

  #rename
  Counts <- cbind(data.s[[1]][genes, ], data.m[[1]][genes, ])
  CountsERCC <- cbind(data.s[[2]], data.m[[2]])
  
  #remove samples from earlier protocol
  keep.plates <- c(
    "NJA01203", "NJA01205","NJA01303", "NJA01401", "NJA01503", "NJA01504",
    "NJA01801", "NJA01803", "NJD00101", "NJD00102", "NJD00103", "NJD00104",
    "NJA01201", "NJA01202", "NJA01301", "NJA01302", "NJA01501", "NJA01901",
    "NJA02001"
  )
  Meta <- filter(Meta, unique_key %in% keep.plates)
  Counts <- Counts[, colnames(Counts) %in% pull(Meta, sample)]
  CountsERCC <- CountsERCC[, colnames(CountsERCC) %in% pull(Meta, sample)]
  
  #save .rda
  if(save) saveRDA(projectName, Counts, CountsERCC, Meta)
}
