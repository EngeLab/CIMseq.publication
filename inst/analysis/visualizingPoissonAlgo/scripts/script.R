packages <- c(
  "CIMseq","CIMseq.data", "CIMseq.testing", "printr", "dplyr",
  "ggthemes", "purrr", "matrixStats", "future", "tibble", "tidyr"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()
print(paste("Running analysis in ", currPath))

#############
#FUNCTIONS
############

syntheticMultipletsFromCounts <- function(
  singlets, classes, fractions,
  seed, idx, ...
){
  set.seed(seed)
  exCounts <- singlets[, idx]
  
  #adjust the counts by multiplying with the corresponding values adjustment
  adjusted <- round(t(t(exCounts) * fractions))
  
  #calculate the sum of counts for each gene
  rs <- matrixStats::rowSums2(adjusted)
  matrix(rs)
}

generateSyntheticMultiplets <- function(
  singlets, classes, fractions,
  seed, n, idx = NULL, ...
){
  mat <- purrr::map(1:n, function(y) {
    syntheticMultipletsFromCounts(
      singlets = singlets, classes = classes,
      fractions = fractions, seed = seed + y, idx = idx[[y]]
    )
  }) %>%
    do.call(cbind, .)
  
  colnames(mat) <- paste0("multiplet.", 1:n)
  rownames(mat) <- rownames(singlets)
  mat
}
############

##RUN ALGO
load('../SCM.analysis/data/CIMseqData.rda')
load('../SCM.analysis/data/sObj.rda')

renameClasses <- function(class) {
  case_when(
    class == "0" ~ "A375",
    class == "1" ~ "HCT116",
    class == "2" ~ "HOS",
    TRUE ~ "error"
  )
}

getData(cObjSng, "classification") <- renameClasses(getData(cObjSng, "classification"))
fractions <- getData(sObj, "fractions")
colnames(fractions) <- renameClasses(colnames(fractions))
sObj@fractions <- fractions

##GENERATE SYNTHETIC MULTIPLETS
#generate synthetic multiplets
g <- c("ANXA3", "CD74", "ACTG2")
class <- getData(cObjSng, "classification")

#A375-HCT116 multiplet: m.NJB00204.A02
m <- "m.NJB00204.A02"
mul <- getData(cObjMul, "counts.cpm")
oneMultiplet <- ceiling(mul[, colnames(mul) == m])
n.multiplets <- 1000
idx <- map(1:n.multiplets, function(n) {
  purrr::map_int(unique(class), ~sample(which(class == .x), 1))
})

#generate synthetic multiplets
#use correct fractions
adjust.right <- as.numeric(getData(sObj, "fractions")[m, unique(class)]) #A375-HCT116
mat <- generateSyntheticMultiplets(
  getData(cObjSng, "counts.cpm")[g, ], class, fractions = adjust.right,
  seed = 8923, n = n.multiplets, idx = idx
) 

synthetic.mul.right <- mat %>%
  matrix_to_tibble("gene") %>%
  add_column(real.multiplet = oneMultiplet[g]) %>%
  gather(sample, value, -gene, -real.multiplet) %>%
  mutate(type = paste0(
    "Correct solution: A375 (", round(adjust.right[3], digits = 2),
    "), HCT116 (", round(adjust.right[2], digits = 2),
    "), HOS (", round(adjust.right[1], digits = 2), ")"
  ))

#use incorrect fractions
adjust.wrong <- rep(1/3, 3)
mat <- generateSyntheticMultiplets(
  getData(cObjSng, "counts.cpm")[g, ], class, fractions = adjust.wrong,
  seed = 8923, n = n.multiplets, idx = idx
)

synthetic.mul.wrong <- mat %>%
  matrix_to_tibble("gene") %>%
  add_column(real.multiplet = oneMultiplet[g]) %>%
  gather(sample, value, -gene, -real.multiplet) %>%
  mutate(type = paste0(
    "Incorrect solution: A375 (", round(adjust.wrong[3], digits = 2),
    "), HCT116 (", round(adjust.wrong[2], digits = 2),
    "), HOS (", round(adjust.wrong[1], digits = 2), ")"
  ))

#bind results into one structure
synthetic.mul <- bind_rows(synthetic.mul.right, synthetic.mul.wrong) %>%
  mutate(gene = case_when(
    gene == "ANXA3" ~ "HCT116: ANXA3",
    gene == "CD74" ~ "A375: CD74",
    gene == "ACTG2" ~ "HOS: ACTG2",
    TRUE ~ "error"
  ))

if(!"data" %in% list.dirs(currPath, full.names = FALSE)) system('mkdir data')
save(synthetic.mul, file = file.path(currPath, "data/syntheticMultiplets.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
