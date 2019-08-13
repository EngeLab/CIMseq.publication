#' palette
#'
#' Figure palette
#'
#' @name palette
#' @rdname palette
#' @param type character; Indicates which order subset to return. Can be "total",
#' "si", or "c".
#' @return A named character vector with names corresponding to classes and 
#' elements corresponding to colours.
#' @author Jason T. Serviss

#' @export
#' @importFrom dplyr case_when

palette <- function(type) {
  out <- case_when(
    type == "total" ~ list(.palette.total()),
    type == "si" ~ list(.palette.si()),
    type == "c" ~ list(.palette.c())
  )
  out[[1]]
}

.palette.total <- function() {
  c(
    `C Stem 1` = "#C2FF99", `C Stem 2` = "#C0B9B2", 
    `C Stem 3` = "#FF8A9A", `C Transit amplifying` = "#6F0062", 
    `C Progenitor` = "#FECFAD",  `C Colonocytes` = "#00489C", 
    `C Goblet proliferating` = "#EEC3FF", `C Goblet Plet1 1` = "#e0a81c", 
    `C Goblet Plet1 2` = "#922329", `C Goblet Plet1-` = "#FFF69F", 
    Tufft = "#ec4339", Enteroendocrine = "#828585", Blood = "#9B9700", 
    `SI Goblet` = "#DDEFFF", `SI Paneth` = "#7B4F4B", `SI Stem` = "#A1C299", 
    `SI Transit amplifying` = "#0AA6D8", `SI Progenitor early` = "#00846F", 
    `SI Progenitor late` = "#FFB500", `SI Enterocytes` = "#C2FFED"
  )
}

.palette.si <- function() {
  c(
    Tufft = "#ec4339", Enteroendocrine = "#828585", Blood = "#9B9700", 
    `SI Goblet` = "#DDEFFF", `SI Paneth` = "#7B4F4B", `SI Stem` = "#A1C299", 
    `SI Transit amplifying` = "#0AA6D8", `SI Progenitor early` = "#00846F", 
    `SI Progenitor late` = "#FFB500", `SI Enterocytes` = "#C2FFED"
  )
}

.palette.c <- function() {
  c(
    `C Stem 1` = "#C2FF99", `C Stem 2` = "#C0B9B2", 
    `C Stem 3` = "#FF8A9A", `C Transit amplifying` = "#6F0062", 
    `C Progenitor` = "#FECFAD",  `C Colonocytes` = "#00489C", 
    `C Goblet proliferating` = "#EEC3FF", `C Goblet Plet1 1` = "#e0a81c", 
    `C Goblet Plet1 2` = "#922329", `C Goblet Plet1-` = "#FFF69F", 
    Tufft = "#ec4339", Enteroendocrine = "#828585", Blood = "#9B9700"
  )
}