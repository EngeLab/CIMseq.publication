
#PACKAGES
packages <- c(
  "CIMseq", "tidyverse", "future", "future.apply", "Seurat"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()
print(paste("Running analysis in ", currPath))

#check package version
algoV <- sessionInfo()$otherPkgs$CIMseq$Version
last3 <- paste(strsplit(algoV, "\\.")[[1]][2:4], collapse = "")
if(!as.numeric(last3) >= 100) {
  stop("CIMseq package version too low. Must be >= 0.2.0.0")
}

currPath <- getwd()

##spCounts
sng.samples <- filter(SCM.Meta, cellNumber == "Singlet" & !filtered)$sample
sng <- SCM.Counts[, sng.samples]
sngERCC <- SCM.CountsERCC[, sng.samples]
mul.samples <- filter(SCM.Meta, cellNumber == "Multiplet" & !filtered)$sample
mul <- SCM.Counts[, mul.samples]
mulERCC <- SCM.CountsERCC[, mul.samples]

#Dimensionality reduction and classification
print(paste0("Starting dim.red and classification analysis at ", Sys.time()))
mca <- CreateSeuratObject(raw.data = sng)
mca <- NormalizeData(
  object = mca, normalization.method = "LogNormalize", scale.factor = 1e6
)
mca <- FindVariableGenes(
  object = mca, mean.function = ExpMean, dispersion.function = LogVMR,
  do.plot = FALSE, x.low.cutoff = 1
)
mca <- ScaleData(
  object = mca, genes.use = mca@var.genes, display.progress = FALSE, do.par = TRUE,
  num.cores = 4
)
mca <- RunPCA(
  object = mca, pc.genes = mca@var.genes, do.print = FALSE
)

mca <- JackStraw(object = mca, num.replicate = 100, display.progress = TRUE, num.pc = 10)
mca <- JackStrawPlot(object = mca, PCs = 1:10)
PCp <- mca@dr$pca@jackstraw@overall.p.values
pcs <- PCp[PCp[, 2] < 0.001, 1]
print(paste0("Using ", max(pcs), " principal components."))

#PCElbowPlot(object = mca, num.pc = 20) + scale_x_continuous(breaks = seq(0, 20, 1))

set.seed(7239832)

mca <- FindClusters(
  object = mca, reduction.type = "pca", dims.use = pcs, resolution = 1,
  save.SNN = TRUE, n.start = 100, nn.eps = 0.5, print.output = FALSE,
  force.recalc = TRUE
)
mca <- RunUMAP(
  object = mca, reduction.use = "pca", dims.use = pcs, min_dist = 0.5,
  n_neighbors = 20
)

# DimPlot(
#   object = mca, reduction.use = "umap", no.legend = FALSE, do.return = TRUE,
#   vector.friendly = FALSE, pt.size = 1
# ) + scale_colour_manual(values = col40())
#
# FeaturePlot(
#   mca, c("CD74", "ANXA3", "ACTG2"),
#   reduction.use = "umap", dark.theme = FALSE, pt.size = 1,
#   vector.friendly = FALSE
# )

markers <- FindAllMarkers(
  object = mca, only.pos = TRUE, min.diff.pct = 0.4, logfc.threshold = log(2),
  test.use = "roc"
)
# n <- min(table(markers$cluster))
# toPlot <- markers %>%
#   group_by(cluster) %>%
#   top_n(max(c(n, 100)), myAUC) %>%
#   pull(gene)
#
# DoHeatmap(
#   object = mca, genes.use = toPlot, slim.col.label = TRUE, remove.key = TRUE,
#   group.label.rot = TRUE
# )
print(paste0("Done all dim.red and classification analysis at ", Sys.time()))

#extract data
if(identical(colnames(sng), rownames(FetchData(mca, "ident")))) {
  classes <- as.character(FetchData(mca, "ident")[[1]])
} else {
  stop("Naming mismatch")
}
select <- match(markers$gene, rownames(sng))
names(select) <- markers$cluster
dim.red <- mca@dr$umap@cell.embeddings
colnames(dim.red) <- NULL

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(sng, sngERCC, dim.red, classes)
cObjMul <- CIMseqMultiplets(mul, mulERCC, select)

#save
if(!"data" %in% list.dirs(currPath, full.names = FALSE)) system('mkdir data')
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))

##spSwarm
baseSeed <- 43892
init <- map(1:10, function(i) {
  cbind(
    swarmInit(cObjSng, 2, null.weight = 0.5, seed = baseSeed + i), 
    swarmInit(cObjSng, 3, null.weight = 0.5, seed = baseSeed + i)
  )
}) %>% do.call(cbind, .)

future::plan(multiprocess)
print(paste0("Starting deconvolution at ", Sys.time()))
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 2, swarmsize = ncol(init), 
  nSyntheticMultiplets = 5, swarmInit = init, 
  psoControl = list(eps.stagnate = 1, maxit.stagnate = 5)
)
print(paste0("Finished deconvolution at ", Sys.time()))

if(!dir.exists('data')) dir.create('data')
save(sObj, file = file.path(currPath, "data/sObj.rda"))

if(!dir.exists('logs')) dir.create('logs')
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))
print("finished")
