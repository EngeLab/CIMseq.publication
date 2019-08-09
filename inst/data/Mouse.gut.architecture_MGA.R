#run from package root
#source('inst/rawData/countsMgfp/Mouse.gut.architecture.R')

packages <- c("dplyr", "stringr")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

MGA <- function(save = TRUE) {
  projectName <- "Mouse.gut.architecture_MGA"
  shortName <- "MGA"
  cat(paste0('Processing ', projectName, '\n'))
  
  m.link <- "https://drive.google.com/open?id=1fl1ZyR5jUJ22dKxfw6Upd5jNmnZqoe5d"
  Meta <- getMetadataViaLink(m.link)
  if("Missing" %in% colnames(Meta)) {
    Meta <- Meta %>%
      filter(is.na(Missing) | Missing == FALSE) %>%
      select(-Missing)
  }
  
  c.link <- c(
    'https://drive.google.com/file/d/1sjFdG8JPR-pvk0IbdJYpCQHB0b8c9jwU/view?usp=sharing',
    'https://drive.google.com/file/d/1om-fGtkJhrSC4n69WgBi2JD8kRrGDXyi/view?usp=sharing',
    'https://drive.google.com/file/d/11rmtUZ9Gyq9jec9fgO3xtkKBZ7LccDmc/view?usp=sharing',
    'https://drive.google.com/file/d/1LldRqRGfRKN8uTDqd0Hk12elTANJS6Kg/view?usp=sharing',
    'https://drive.google.com/file/d/1bSGsxNdajABpUg_wDUHL6eNLJjwpP9Us/view?usp=sharing',
    'https://drive.google.com/file/d/131WFJmLFMZlnnPQjMPijyMuXulo2ZWoN/view?usp=sharing',
    'https://drive.google.com/file/d/1R8d1AEkOnHR4GCsbDlvNF10u1qqHA0tz/view?usp=sharing',
    'https://drive.google.com/file/d/1J6VHAGQVZOFT6pzD6oMtzloqTKQl94jl/view?usp=sharing',
    'https://drive.google.com/file/d/1SQnbGLcee3E4m7uT-yqVzk0DrJtV2zz_/view?usp=sharing',
    'https://drive.google.com/file/d/1G-X3xuGsvWGm5lA0HskFN8rNdvkjj_Jp/view?usp=sharing',
    'https://drive.google.com/file/d/1Tzrw7kHaLLfaPzjoQyWf5sZT97LVsVKP/view?usp=sharing',
    'https://drive.google.com/file/d/1HIE-vTB1NInrACfbIKJe8iMX6LHyfmxm/view?usp=sharing',
    'https://drive.google.com/file/d/1VFZFJLtSQjXlzCp-JsrG6S3t28vDxFHf/view?usp=sharing',
    'https://drive.google.com/file/d/1evrm5lGjfOoNY6mb63Igmp5pljzyJPv-/view?usp=sharing'
  )
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
