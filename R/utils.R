#' seuratToCIMseq
#'
#' Exttracts info from Seurat object and returns CIMseqSinglets and 
#' CIMseqMultiplets objects.
#'
#' @name seuratToCIMseq
#' @rdname seuratToCIMseq
#' @param seurat.obj Seurat; A Seurat object.
#' @param singlets matrix; A matrix with singlet counts.
#' @param singlet.ercc matrix; A matrix with singlet ERCC counts.
#' @param multiplets matrix; A matrix with multiplet counts.
#' @param multiplets matrix; A matrix with multiplet ERCC counts.
#' @param markers data.frame; Output from \code{\link[Seurat]{FindAllMarkers}}
#'  function. 
#' @return A list with CIMseqSinglets object as the first element and 
#' CIMseqMultiplets object as the second element.
#' @author Jason T. Serviss
#'
NULL
#' @export
#' @import CIMseq

seuratToCIMseq <- function(
  seurat.obj, singlets, singlet.ercc, multiplets, multiplet.ercc, markers
){
  #singlet data
  singlets <- singlets[, colnames(singlets) %in% colnames(mca@data)]
  singletERCC <- singletERCC[, colnames(singletERCC) %in% colnames(singlets)]
  idx <- match(names(mca@ident), colnames(singlets))
  
  #classifications
  classes <- as.character(mca@ident)[idx]
  names(classes) <- names(mca@ident)[idx]
  
  #features
  var.genes <- unique(markers$gene)
  select <- which(rownames(singlets) %in% var.genes)
  
  #dim red
  dim.red <- mca@dr$umap@cell.embeddings
  colnames(dim.red) <- NULL
  
  #CIMseq objs
  cObjSng <- CIMseqSinglets(singlets, singletERCC, dim.red, classes)
  cObjMul <- CIMseqMultiplets(multiplets, multipletERCC, select)
  
  return(list(cObjSng, cObjMul))
}