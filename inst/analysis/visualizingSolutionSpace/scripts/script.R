packages <- c(
  "CIMseq", "printr", "ggthemes", "tibble", "future",
  "purrr", "dplyr", "tidyr"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()
print(paste("Running analysis in ", currPath))

load('../SCM.analysis/data/CIMseqData.rda')
keep <- c("m.NJB00204.G04", "m.NJB00204.D07", "m.NJB00204.F12")
cObjMul <- CIMseqMultiplets(
  getData(cObjMul, "counts")[, keep],
  getData(cObjMul, "counts.ercc")[, keep],
  getData(cObjMul, "features")
)

plan(multiprocess)
sObj <- CIMseqSwarm(
  cObjSng, cObjMul, maxiter = 40, swarmsize = 200,
  nSyntheticMultiplets = 500, report = TRUE, reportRate = 1
)

ss <- getData(sObj, "fractions")
colnames(ss) <- c("HOS", "HCT116", "A375")
sObj@fractions <- ss

if(!"data" %in% list.dirs(currPath, full.names = FALSE)) system('mkdir data')
save(cObjSng, cObjMul, sObj, file = file.path(currPath, 'data/algoOutput.rda'))

renameClasses <- function(class) {
  case_when(
    class == "0" ~ "A375", class == "1" ~ "HCT116", class == "2" ~ "HOS",
    TRUE ~ "error"
  )
}
report <- unnest(getData(sObj, "stats"))
colnames(report)[6:8] <- renameClasses(colnames(report)[6:8])

save(report, file = file.path(currPath, "data/report.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo.txt"))

