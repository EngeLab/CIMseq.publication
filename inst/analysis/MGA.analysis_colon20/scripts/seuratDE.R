library("Seurat")

currPath <- getwd()

#load data
path1 <- 'data/seurat.rda'

if(file.exists(path1)) {
  load(path1)
}
rm(list=c("path1"))

#run DE
seurat <- SubsetData(mca.c, ident.use = c("1", "8", "10"))
seurat <- RenameIdent(seurat, "1", "Goblet Plet1 neg.")
seurat <- RenameIdent(seurat, "8", "Goblet Plet1 pos.")
seurat <- RenameIdent(seurat, "10", "Goblet Plet1 pos.")
seuratDE <- FindMarkers(seurat, ident.1 = "Goblet Plet1 pos.", only.pos = TRUE)

#Save
save(seuratDE, file = file.path(currPath, "data/seuratDE.rda"))
