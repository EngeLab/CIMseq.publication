library(purrr)
library(CIMseq.publication)

#ARGUMENTS
directories <- c(
  "MGA.analysis_enge20",
  "SCM.analysis", 
  "deconvoluteSingletsMouse", 
  "MGA.analysis_colon20", 
  "MGA.analysis_SI20", 
  "multipletCellNrDist", 
  "visualizingPoissonAlgo", 
  "visualizingSolutionSpace"
)
directories <- file.path('./inst/analysis', directories)


allFigures <- list(
  "Figure 1" = c(
    b = 'inst/analysis/SCM.analysis/figures/SCM.ercc.pdf', 
    c = 'inst/analysis/SCM.analysis/figures/SCM.estimatedVSactual.pdf',
    d = 'inst/analysis/multipletCellNrDist/figures/multipletCellNrDist.pdf',
    e = 'inst/analysis/SCM.analysis/figures/figure2.pdf'
  ), 
  "Figure 2" = c(
    a1 = 'inst/analysis/MGA.analysis_SI20/figures/MGA.enge20SI.markers.pdf',
    a2 = 'inst/analysis/MGA.analysis_SI20/figures/MGA.enge20SI.classes.pdf',
    a3 = 'inst/analysis/MGA.analysis_SI20/figures/MGA.enge20SI.cellcycle.pdf',
    b = 'inst/analysis/MGA.analysis_SI20/figures/MGA.enge20SI.sngMulRelFreq.pdf',
    c = 'inst/analysis/MGA.analysis_SI20/figures/MGA.enge20SI.circos.pdf'
  ), 
  "Figure 3" = c(
    a1 = 'inst/analysis/MGA.analysis_colon20/figures/MGA.enge20.markers.pdf',
    a2 = 'inst/analysis/MGA.analysis_colon20/figures/MGA.enge20.classes.pdf',
    a3 = 'inst/analysis/MGA.analysis_colon20/figures/MGA.enge20.Mki67.pdf',
    b = 'inst/analysis/MGA.analysis_colon20/figures/MGA.enge20.sngMulRelFreq.pdf',
    c = 'inst/analysis/MGA.analysis_colon20/figures/MGA.enge20.circos.pdf'
  ), 
  "Figure 4" = c(
    a = 'inst/analysis/MGA.analysis_enge20/figures/MGA.enge.stemness.dotplot.pdf'
  )
)

#FUNCTIONS
#generate data
generateData <- function(directory) {
  home <- getwd()
  print(paste0("Processing ", basename(directory)))
  setwd(directory)
  scripts <- list.files('scripts', full.names = TRUE)
  if(length(scripts) > 0) {
    for(y in 1:length(scripts)){
      print(paste0("Running ", scripts[y]))
      source(scripts[y])
    }
  }
  setwd(home)
}

#re-run analysis
runAnalysis <- function(directory) {
  home <- getwd()
  print(paste0("Processing ", basename(directory)))
  setwd(directory)
  analysis <- list.files('analysis', full.names = TRUE, pattern = ".Rmd")
  for(j in 1:length(analysis)) {
    rmarkdown::render(analysis[j])
  }
  setwd(home)
}

#generate figures
generateFigures <- function(figureSpecs, out, name) {
  first <- "```{r, fig.align='center', out.width='450px', out.height='450px'}\n"
  last = "```\n"
  head <- paste0(
    '---\ntitle: ', name, '\nauthor: "Jason T. Serviss"\noutput:\n  ',
    'html_document:\n    toc: yes\n    code_folding: hide\n---\n\n'
  )
  
  tmp <- tempfile(fileext = ".Rmd")
  cat(head, file = tmp)
  for(f in figureSpecs) {
    cat(first, file = tmp, append = TRUE)
    cat(paste0("knitr::include_graphics('", f, "')\n"), file = tmp, append = TRUE)
    cat(last, file = tmp, append = TRUE)
    cat('<br></br>\n', file = tmp, append = TRUE)
  }
  
  rmarkdown::render(
    tmp, output_format = "html_document", output_dir = out, output_file = paste0(name, ".html")
  )
}

#RUN
one <- map(directories, ~generateData(.x))
two <- map(directories, ~runAnalysis(.x))
three <- map(1:length(allFigures), ~generateFigures(
  file.path(getwd(), allFigures[[.x]]),
  'inst/figures',
  names(allFigures)[.x])
)