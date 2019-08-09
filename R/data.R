#' getCountsDataViaLinks
#'
#' Downloads the count data project and concatenates it into a data.frame. 
#' The column containing the gene name, to be named HGN, is expected in all 
#' counts files.
#'
#' @name getCountsDataViaLinks
#' @rdname getCountsDataViaLinks
#' @param links character; The links to the counts files to download.
#' @return A data.frame with the unfiltered counts data.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom googledrive drive_get drive_download drive_auth_config as_id
#' @importFrom purrr map2 map reduce
#' @importFrom dplyr full_join
#' @importFrom utils read.table

getCountsDataViaLinks <- function(links) {
  drive_auth_config(active = FALSE)
  countsFiles <- drive_get(as_id(links))$name
  localPaths <- file.path(tempdir(), countsFiles)
  
  #download
  trash <- purrr::map2(links, localPaths, function(fp, lp) {
    googledrive::drive_download(as_id(fp), lp, overwrite = TRUE)
  })
  rm(trash)
  
  #load counts
  loaded <- purrr::map(localPaths, utils::read.table, header = TRUE, sep = "\t")
  purrr::reduce(loaded, dplyr::full_join, by = "HGN")
}

#' getMetadataViaLink
#'
#' Utilizes the EngeMetadata package to download and format the metadata
#' associated with a specific project folder. Downloads all metadata/plate info
#' present in the annotation folder by default.
#'
#' @name getMetadataViaLink
#' @rdname getMetadataViaLink
#' @param link character; Link to the processed metadata.
#' @return The metadata tibble containing the concatenated information in the 
#' project annotation folder.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom googledrive drive_get drive_download as_id
#' @importFrom readr read_tsv
#' @importFrom tibble as_tibble

getMetadataViaLink <- function(link) {
  drive_auth_config(active = FALSE)
  
  name <- googledrive::drive_get(as_id(link))$name
  localpath <- file.path(tempdir(), name)
  googledrive::drive_download(as_id(link), localpath, overwrite = TRUE)
  as_tibble(read.table(localpath, stringsAsFactors = FALSE))
}

#' Rename GFP mouse samples.
#'
#' Updates sample names from the old nomenclature (GFPpos.C1.Doublet.5.F9.htseq)
#' to the new nomenclature (NJA00103.D09.htseq).
#'
#' @name renameMgfpSamples
#' @rdname renameMgfpSamples
#' @param oldNames character; A character vector of the old names that have
#'  already been processed by the removeHTSEQsuffix functions. Names that are
#' already in the new nomenclature will not be effected.
#' @return A character vector of the updated names.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom stringr str_detect str_replace str_extract
#' @importFrom dplyr case_when

renameMgfpSamples <- function(oldNames) {
  bool1 <- str_detect(oldNames, "GFPneg.C1.Singlet.4")
  bool2 <- str_detect(oldNames, "GFPneg.SI1.Singlet.2")
  bool3 <- str_detect(oldNames, "GFPpos.C1.Doublet.5")
  bool4 <- str_detect(oldNames, "GFPpos.C1.Singlet.4")
  bool5 <- str_detect(oldNames, "GFPpos.SI1.Doublet.3")
  bool <- bool1 | bool2 | bool3 | bool4 | bool5
  idx <- which(bool)
  
  #extract and fix plate positions
  pos <- paste0(".", str_extract(oldNames[bool], "[^.]+$"))
  posidx <- which(nchar(pos) == 3)
  pos[posidx] <- str_replace(pos[posidx], "(..)(.)", "\\10\\2")
  
  #update plate info
  plate <- case_when(
    bool1 ~ "NJA00110",
    bool2 ~ "NJA00101",
    bool3 ~ "NJA00111",
    bool4 ~ "NJA00110",
    bool5 ~ "NJA00107",
    TRUE ~ oldNames
  )
  
  #synthesize new names
  plate[bool] <- paste0(plate[bool], pos)
  plate
}

#' moveGenesToRownames
#'
#' Moves the first column of the counts data.frame to rownames and removes
#' the old column.
#'
#' @name moveGenesToRownames
#' @rdname moveGenesToRownames
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame is returned with the first column moved to
#' the rownames of the data.frame.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(LETTERS, a = runif(26, 1, 100))
#' moveGenesToRownames(counts)
#'
NULL
#' @export

moveGenesToRownames <- function(counts) {
  if(!length(dim(counts))) {
    stop("dim(counts) must have a positive length.")
  }
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if((class(counts[[1]]) != "character") & (class(counts[[1]]) != "factor")) {
    stop("The first column of counts is not class character or class factor.")
  }
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  return(counts)
}

#' removeHTSEQsuffix
#'
#' HTSeq adds the suffix ".htseq" to column names when it reports counts. This
#' function removes that suffix from the column names of the supplied counts
#' data.frame.
#'
#' @name removeHTSEQsuffix
#' @rdname removeHTSEQsuffix
#' @param data data.frame or character; A data frame with sample IDs as colnames
#'  or a vector of sample IDs.
#' @return The counts data.frame with the ".htseq" suffix removed from the
#'  column names.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a.htseq = runif(26), b.htseq = runif(26))
#' removeHTSEQsuffix(counts)
#'
NULL
#' @export

removeHTSEQsuffix <- function(data) {
  .inputChecks_removeHTSEQsuffix(data)
  if(class(data) == "data.frame") {
    .dataframeChecks_removeHTSEQsuffix(data)
    colnames(data) <- gsub("(.*)\\.htseq$", "\\1", colnames(data))
    return(data)
  }
  if(class(data) == "character") {
    data <- gsub("(.*)\\.htseq$", "\\1", data)
    return(data)
  }
}

.inputChecks_removeHTSEQsuffix <- function(data) {
  if(!class(data) %in% c("data.frame", "character")) {
    stop("The data argument must be a character vector or data.frame.")
  }
}

.dataframeChecks_removeHTSEQsuffix <- function(data) {
  if(is.null(colnames(data))) {
    stop("is.null(colnames(data)) returned TRUE.")
  }
  if(all(colnames(data) == paste0("V", 1:ncol(data)))) {
    warning("Your colnames are V1..Vi. Are these sample names?" )
  }
  if(!any(grepl("\\.htseq", colnames(data)))) {
    warning("Could not find the .htseq suffix in colnames(data)")
  }
}

#' labelSingletsAndMultiplets
#'
#' Adds the prefix "s." to colnames of columns that contain singlets and "m." to
#' columns that contain multiplets.
#'
#' @name labelSingletsAndMultiplets
#' @rdname labelSingletsAndMultiplets
#' @param data data.frame or character; A data frame with sample IDs as colnames
#'  or a vector of sample IDs.
#' @param ids character; a character vector of regex statments used to indicate
#'  colnames of columns that contain SINGLETS.
#' @return If data is of class data.frame colnames are modified in place with
#' the appropriate prefix attatched. If data is of class character, the modified
#' character vector is returned. Throws a warning if no singlets are detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26), b = runif(26), c = runif(26))
#' labelSingletsAndMultiplets(counts, ids = c("a", "b"))
#'
NULL
#' @export

labelSingletsAndMultiplets <- function(data, ids) {
  .inputChecks_labelSingletsAndMultiplets(data, ids)
  
  if(class(data) == "data.frame") {
    
    .dataframeChecks_labelSingletsAndMultiplets(data)
    bool <- sapply(ids, function(x) grepl(x, colnames(data)))
    .boolChecks_labelSingletsAndMultiplets(bool)
    colnames(data) <- ifelse(
      rowSums(bool) > 0,
      paste("s.", colnames(data), sep = ""),
      paste("m.", colnames(data), sep = "")
    )
    return(data)
    
  } else {
    
    bool <- sapply(ids, function(x) grepl(x, data))
    .boolChecks_labelSingletsAndMultiplets(bool)
    data <- ifelse(
      rowSums(bool) > 0,
      paste("s.", data, sep = ""),
      paste("m.", data, sep = "")
    )
    return(data)
  }
}

.inputChecks_labelSingletsAndMultiplets <- function(data, ids) {
  if(class(ids) != "character") {
    stop("The id argument must be a character vector.")
  }
  if(!class(data) %in% c("data.frame", "character")) {
    stop("The data argument must be a character vector or data.frame.")
  }
  if(any(is.na(ids))) {
    stop("The ids argument contains NA values.")
  }
}

.dataframeChecks_labelSingletsAndMultiplets <- function(data) {
  if(is.null(colnames(data))) {
    stop("is.null(colnames(data)) returned TRUE.")
  }
  if(!all(sapply(data, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
}

.boolChecks_labelSingletsAndMultiplets <- function(bool) {
  if(sum(bool) == 0) {
    warning("No singlets were detected. Is your id argument valid?")
  }
}

#' detectERCCreads
#'
#' Detects the 92 ERCC reads in the rownames of the counts data.frame. ERCC
#' reads must be named with the standard naming convention that matches the
#' regex "^ERCC\\-[0-9]*$".
#'
#' @name detectERCCreads
#' @rdname detectERCCreads
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains an ERCC read. A warning is issued if all 92
#' ERCC reads are not detected.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(100), row.names = c(1:8, paste0("ERCC-", 1:92)))
#' detectERCCreads(counts)
#'
NULL
#' @export

detectERCCreads <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }
  ercc <- grepl("^ERCC\\-[0-9]*$", rownames(counts))
  
  if(sum(ercc) != 92) {
    warning("Couldn't detect all ERCC reads.")
  }
  
  return(ercc)
}

#' detectNonGenes
#'
#' HTSeq outputs several "non-genes" in the counts output. These include:
#' "__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned",
#' "__alignment_not_unique". Since these are measurements of the quantification
#' quality and are often undesirable in downstream analysis, this function
#' detects them in the submitted counts data.frame to allow easy removal.
#'
#' @name detectNonGenes
#' @rdname detectNonGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row contains a non-gene.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(runif(10), row.names = c(1:9, "__no_feature"))
#' detectNonGenes(counts)
#'
NULL
#' @export

detectNonGenes <- function(counts) {
  if(class(rownames(counts)) != "character") {
    stop("rownames(counts) is not of class character.")
  }
  if(all(rownames(counts) == as.character(1:nrow(counts)))) {
    m <- "rownames(counts) = 1:nrow(counts). Are gene names in rownames counts?"
    stop(m)
  }
  
  nonGenes <- c(
    "__no_feature", "__ambiguous", "__too_low_aQual",
    "__not_aligned", "__alignment_not_unique"
  )
  rownames(counts) %in% nonGenes
}

#' filterCountsData
#'
#' This is a wrapper function for \code{\link{detectLowQualityGenes}},
#' \code{\link{detectLowQualityCells.totalCounts}},
#' \code{\link{detectLowQualityCells.ERCCfrac}},
#' \code{\link{detectLowQualityCells.housekeeping}},
#' \code{\link{convertCountsToMatrix}}. It uses the
#' "cellNumber" column in the metadata to annotate singlets and multiplets and,
#' therefore, the column must be present in the metadata. In addition, singlets
#' should be named "Singlet" in the metadata. The function also checks for the
#' presence of NAs in the counts data and returns an error if NA values are
#' detected.
#'
#' @name filterCountsData
#' @rdname filterCountsData
#' @param counts data.frame; The formated counts data.
#' @param countsERCC data.frame; The formated counts ERCC data.
#' @param filters Character; A vector indicating the types of filtering to be
#' performed. Options include: "genes", "totalCounts",
#' "ERCCfrac", and "housekeeping".
#' @param geneMinCount See \code{\link{detectLowQualityGenes}}
#' @param cellMinCount See \code{\link{detectLowQualityCells.totalCounts}}
#' @param geneName See \code{\link{detectLowQualityCells.housekeeping}}
#' @param quantileCut See \code{\link{detectLowQualityCells.housekeeping}}
#' @param percentile See \code{\link{detectLowQualityCells.ERCCfrac}}.
#' @return A list with the filtered counts as the first element and the filtered
#' ERCC genes as the second element.
#' @author Jason Serviss
#'
NULL
#' @export

filterCountsData <- function(
  counts, countsERCC,
  filters = c("genes", "totalCounts", "ERCCfrac", "housekeeping"),
  geneMinCount, cellMinCount, geneName, quantileCut.hk, quantileCut.ercc
){
  #remove low quality genes
  if("genes" %in% filters) {
    counts <- counts[detectLowQualityGenes(counts, mincount = geneMinCount), ]
  }
  
  #remove low quality cells
  lqc.totalCounts <- lqc.housekeeping <- lqc.ERCCfrac <- rep(TRUE, ncol(counts))
  
  if("totalCounts" %in% filters) {
    lqc.totalCounts <- detectLowQualityCells.totalCounts(
      counts,
      mincount = cellMinCount
    )
  }
  
  if("housekeeping" %in% filters) {
    lqc.housekeeping <- detectLowQualityCells.housekeeping(
      counts,
      geneName = geneName,
      quantileCut = quantileCut.hk
    )
  }
  
  if("ERCCfrac" %in% filters) {
    lqc.ERCCfrac <- detectLowQualityCells.ERCCfrac(
      counts,
      countsERCC,
      quantileCut = quantileCut.ercc
    )
  }
  
  lqc <- lqc.totalCounts & lqc.housekeeping & lqc.ERCCfrac
  
  print(paste0(
    "Removing a total of ",
    sum(!lqc),
    " cells based on the calculated metrics."
  ))
  
  counts <- counts[, lqc]
  countsERCC <- countsERCC[, lqc]
  
  #coerce to matrix
  counts <- convertCountsToMatrix(counts)
  countsERCC <- convertCountsToMatrix(countsERCC)
  
  return(list(counts, countsERCC))
}

#' detectLowQualityGenes
#'
#' In gene expression counts data if can often be the case that some genes are
#' not detected. This can simply be due to the fact that the gene is not
#' expressed in the tissue or, in addition, that the sequencing depth was not
#' sufficient to detect the gene. In addition, some genes may be detected but in
#' so few samples or at such a low level that it makes the quantified value
#' highly unreliable. In these cases, it is desireable to remove the gene before
#' downstream analysis which is facilitated by this function.
#'
#' @name detectLowQualityGenes
#' @rdname detectLowQualityGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @param mincount numeric; A minimum rowSum for which rows with a higher rowSum
#' will be detected. Default = 0.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row meets both tested conditions.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(c(0, runif(10)), c(0, runif(10)), c(0, runif(10)))
#' detectLowQualityGenes(counts)
#'
NULL
#' @export

detectLowQualityGenes <- function(
  counts,
  mincount = 0
){
  #input checks
  
  bool <- rowSums(counts) > mincount
  
  message <- paste0(
    "Detected ", sum(!bool), " low quality genes out of ", nrow(counts),
    " genes input (", round(100 * (sum(!bool) / nrow(counts)), digits = 2),
    "%)."
  )
  print(message)
  return(bool)
}

#' detectLowQualityCells.totalCounts
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by detecting cells
#' whose sum across all genes is > mincount.
#'
#' @name detectLowQualityCells.totalCounts
#' @rdname detectLowQualityCells.totalCounts
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param mincount numeric; A minimum colSum for which columns with a higher
#' colSum will be detected. Default = 4e5.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with colSums > mincount.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detectLowQualityCells.totalCounts(counts)
#'
NULL
#' @export
#' @importFrom stats quantile

detectLowQualityCells.totalCounts <- function(
  counts,
  mincount = 4e5
){
  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)
  
  #colsums check
  cs <- colSums(counts) > mincount
  output[cs] <- TRUE
  
  if(sum(cs) < 2) {
    stop("One or less samples passed the colSums check.")
  }
  
  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on total counts."
  )
  print(message)
  return(output)
}

#' detectLowQualityCells.housekeeping
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by utilizing a house
#' keeping gene and assuming its log2 expression to be normally distributed.
#' We then detect samples where the probability of the expression for the house
#' keeping gene in that sample is greater than the quantile.cut argument.
#'
#' @name detectLowQualityCells.housekeeping
#' @rdname detectLowQualityCells.housekeeping
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param geneName character; The gene name to use for the quantile cutoff. This
#' must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantileCut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.01.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with meeting the criteria specified
#' by the arguments.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1.5, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detectLowQualityCells.housekeeping(counts, geneName = "a1")
#'
NULL
#' @export
#' @importFrom stats median qnorm

detectLowQualityCells.housekeeping <- function(
  counts,
  geneName = 'ACTB',
  quantileCut = 0.01
){
  #input checks
  ##check that geneName is in rownames counts
  if(!geneName %in% rownames(counts)) {
    stop("geneName is not found in rownames(counts)")
  }
  
  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)
  
  #house keeping check
  counts.log <- .norm.log.counts(counts)
  cl.act <- counts.log[geneName, colSums(counts) != 0]
  cl.act.m <- median(cl.act)
  cl.act.sd <- sqrt(
    sum((cl.act[cl.act > cl.act.m] - cl.act.m) ^ 2) /
      (sum(cl.act > cl.act.m) - 1)
  )
  my.cut <- qnorm(p = quantileCut, mean = cl.act.m, sd = cl.act.sd)
  bool <- counts.log[geneName, ] >= my.cut
  output[bool] <- TRUE
  
  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on ", geneName, " expression."
  )
  print(message)
  return(output)
}

#calculates log cpm
.norm.log.counts <- function(counts) {
  log2(.norm.counts(counts) + 1)
}

#calculates cpm
.norm.counts <- function(counts) {
  t(t(counts) / colSums(counts) * 10^6)
}

#' detectLowQualityCells.ERCCfrac
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by indicating samples
#' that have 0 ERCC reads and those that have a fraction of ERCC reads .
#'
#' @name detectLowQualityCells.ERCCfrac
#' @rdname detectLowQualityCells.ERCCfrac
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param geneName character; The gene name to use for the quantile cutoff. This
#' must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantileCut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.99.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with meeting the criteria specified
#' by the arguments.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1.5, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detectLowQualityCells.housekeeping(counts, geneName = "a1")
#'
NULL
#' @export

detectLowQualityCells.ERCCfrac <- function(
  counts,
  ercc,
  quantileCut = 0.99
){
  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)
  
  #calculate fraction of ercc
  cs <- colSums(counts)
  cs.ercc <- colSums(ercc)
  frac.ercc <-  cs.ercc / (cs.ercc + cs)
  frac.ercc[is.na(frac.ercc)] <- 0 #when cs == 0
  l.frac.ercc <- log2(frac.ercc)
  
  #calculate normal approximation and cut
  #median less sensitive to outliers
  #is.infinite needed for when cs.ercc == 0
  m <- median(l.frac.ercc[!is.infinite(l.frac.ercc)], na.rm = TRUE)
  sd <- sd(l.frac.ercc[!is.infinite(l.frac.ercc)], na.rm = TRUE)
  p.cut.h <- 2^qnorm(quantileCut, m, sd)
  p.cut.l <- 2^qnorm(1 - quantileCut, m, sd)
  bool.h <- frac.ercc <= p.cut.h
  bool.l <- frac.ercc >= p.cut.l
  output[bool.h & bool.l] <- TRUE
  
  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on ERCC fraction."
  )
  print(message)
  return(output)
}

#' convertCountsToMatrix
#'
#' Coerces the counts data.frame into a matrix.
#'
#' @name convertCountsToMatrix
#' @rdname convertCountsToMatrix
#' @param counts data.frame; A data frame with counts data.
#' @return The counts data.frame coerced into a matrix.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(a = runif(26, 1, 100), b = runif(26, 1, 100))
#' convertCountsToMatrix(counts)
#'
NULL
#' @export

convertCountsToMatrix <- function(counts) {
  if(class(counts) != "data.frame") {
    stop("Counts is not a data.frame")
  }
  if(!all(sapply(counts, class) %in% c("numeric", "integer"))) {
    stop("Non-numeric columns detected. Should gene names be moved to rownames?")
  }
  as.matrix(counts)
}

#' saveRDA
#'
#' Renames and saves filtered/formated project data as a .rda file in the "data"
#' folder of the current directory.
#'
#' @name saveRDA
#' @rdname saveRDA
#' @param projectName character; The name of the project. Must correspond with
#' the project folder name located at EngeLab/data.
#' @param ... Items to be saved in the .rda file. Note that these variables are
#' renamed before saving such that they have the projectName argument prepended
#' to the current variable name.
#' @return The .rda file is saved as "data/projectName.rda".
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom stringr str_replace

saveRDA <- function(projectName, ...) {
  data <- list(...)
  shortName <- str_replace(projectName, ".*_(.*)", "\\1")
  names <- sapply(substitute(list(...))[-1], deparse)
  names(data) <- paste(shortName, names, sep = ".")
  
  #reassign the variable name from "counts" etc. to "projectName Counts" etc.
  for(i in 1:length(data)) {
    assign(names(data)[i], data[[i]])
  }
  
  save(
    list = names(data),
    file = file.path("./data", paste0(projectName, ".rda")),
    compress = "bzip2"
  )
}
