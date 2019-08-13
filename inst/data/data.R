#source('./inst/data/data.R')

#run scripts to generate data/.rda files
processRaw <- function(
  basePath = './inst/data', ignore = './inst/data/data.R', save = TRUE
){
  print("Downloading and processing data")
  paths <- list.files(
    path = basePath, recursive = TRUE, pattern = '\\.R$', full.names = TRUE
  )

  if(!is.null(ignore)) paths <- paths[!paths %in% ignore]
  funNames <- gsub(".*_(.*)\\.R", "\\1", basename(paths))

  trash <- purrr::map2(paths, funNames, function(p, fn) {
    source(p)
    c.fn <- get(fn)
    c.fn(save = save)
    gc()
  })
}

processRaw()