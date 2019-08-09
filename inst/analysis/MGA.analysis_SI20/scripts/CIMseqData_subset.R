packages <- c("CIMseq", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

#load data
path1 <- '../MGA.analysis_enge20/data/seurat.rda'
path2 <- '../MGA.analysis_enge20/data/CIMseqData.rda'

if(file.exists(path1)) {
  load(path1)
}
if(file.exists(path2)) {
  load(path2)
}
rm(list=c("path1", "path2"))

#subset small intestine samples
keep.plates.SI <- c(
  "NJA01201", "NJA01202", "NJA01301", "NJA01302", "NJA01501"
)
sng.samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates.SI & cellNumber == "Singlet")$sample
stopifnot(all(sng.samples %in% colnames(mca@raw.data)))
mca.si <- SubsetData(mca, sng.samples)

#detect features in small intestine subset
markers <- FindAllMarkers(
  object = mca.si, only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = log(2),
  test.use = "roc"
)
select <- which(rownames(getData(cObjSng, "counts")) %in% unique(markers$gene))

#setup CIMseqData objects
cObjSng <- CIMseqSinglets(
  getData(cObjSng, "counts")[, sng.samples],
  getData(cObjSng, "counts.ercc")[, sng.samples],
  getData(cObjSng, "dim.red")[sng.samples, ],
  getData(cObjSng, "classification")[sng.samples]
)

mul.samples <- filter(MGA.Meta, !filtered & unique_key %in% keep.plates.SI & cellNumber == "Multiplet")$sample
cObjMul <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, mul.samples],
  getData(cObjMul, "counts.ercc")[, mul.samples],
  select
)

#save
save(cObjSng, cObjMul, file = file.path(currPath, "data/CIMseqData.rda"))
save(mca.si, markers, file = file.path(currPath, "data/seurat.rda"))

#write logs
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_CIMseqData.txt"))