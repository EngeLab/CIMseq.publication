#' renameClasses.MGA
#'
#' MGA dataset class renaming (from Seurat default) including tissue 
#' abbreviation.
#'
#' @name renameClasses.MGA
#' @rdname renameClasses.MGA
#' @param class character; Old class names.
#' @param abrev logical; Should tissue abbreviation be included?
#' @return A character vector with new class names.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr case_when

renameClasses.MGA <- function(class, abrev = TRUE) {
  if(abrev) {
    case_when(
      class == "0" ~ "C Stem 1",
      class == "1" ~ "C Goblet Plet1 neg.",
      class == "2" ~ "C Colonocytes",
      class == "3" ~ "SI Stem",
      class == "4" ~ "SI Transit amplifying",
      class == "5" ~ "C Stem 3",
      class == "6" ~ "C Progenitor",
      class == "7" ~ "C Transit amplifying",
      class == "8" ~ "C Goblet Plet1 1",
      class == "9" ~ "SI Goblet",
      class == "10" ~ "C Goblet Plet1 2",
      class == "11" ~ "SI Progenitor late",
      class == "12" ~ "SI Progenitor early",
      class == "13" ~ "C Stem 2",
      class == "14" ~ "Enteroendocrine",
      class == "15" ~ "Tufft",
      class == "16" ~ "SI Enterocytes",
      class == "17" ~ "SI Paneth",
      class == "18" ~ "C Goblet proliferating",
      class == "19" ~ "Blood",
      TRUE ~ "error"
    )
  } else {
    case_when(
      class == "0" ~ "Stem 1",
      class == "1" ~ "Goblet Plet1 neg.",
      class == "2" ~ "Colonocytes",
      class == "3" ~ "Stem",
      class == "4" ~ "Transit amplifying",
      class == "5" ~ "Stem 3",
      class == "6" ~ "Progenitor",
      class == "7" ~ "Transit amplifying",
      class == "8" ~ "Goblet Plet1 1",
      class == "9" ~ "Goblet",
      class == "10" ~ "Goblet Plet1 2",
      class == "11" ~ "Progenitor late",
      class == "12" ~ "Progenitor early",
      class == "13" ~ "Stem 2",
      class == "14" ~ "Enteroendocrine",
      class == "15" ~ "Tufft",
      class == "16" ~ "Enterocytes",
      class == "17" ~ "Paneth",
      class == "18" ~ "Goblet proliferating",
      class == "19" ~ "Blood",
      TRUE ~ "error"
    )
  }
}

#' classOrder.MGA
#'
#' Returns the desired class order for the MGA dataset.
#'
#' @name classOrder.MGA
#' @rdname classOrder.MGA
#' @param type character; Indicates which order subset to return. Can be "total",
#' "si", or "c".
#' @return A character vector with class order
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr case_when

classOrder.MGA <- function(type) {
  out <- case_when(
    type == "total" ~ list(.classOrder.MGA.total()),
    type == "si" ~ list(.classOrder.MGA.SI()),
    type == "c" ~ list(.classOrder.MGA.C())
  )
  out[[1]]
}

.classOrder.MGA.total <- function() {
  c(
    "SI Paneth", "SI Stem", "SI Transit amplifying", "SI Progenitor early", 
    "SI Progenitor late", "SI Enterocytes", "SI Goblet", 
    "C Goblet proliferating", "C Goblet Plet1 1", "C Goblet Plet1 2", "C Stem 1",
    "C Stem 2", "C Stem 3", "C Transit amplifying", "C Progenitor", 
    "C Colonocytes", "C Goblet Plet1 neg.", "Enteroendocrine", "Tufft", "Blood"
  )
}

.classOrder.MGA.SI <- function() {
  c(
    "Paneth", "Stem", "Transit amplifying", "Progenitor early", 
    "Progenitor late", "Enterocytes", "Goblet", "Enteroendocrine", 
    "Tufft", "Blood"
  )
}

.classOrder.MGA.C <- function() {
  c(
    "Goblet proliferating", "Goblet Plet1 1", "Goblet Plet1 2", 
    "Stem 1", "Stem 2", "Stem 3", "Transit amplifying", 
    "Progenitor", "Colonocytes", "Goblet Plet1 neg.",  "Enteroendocrine", 
    "Tufft", "Blood"
  )
}

#' renameClasses.SCM
#'
#' SCM dataset class renaming (from Seurat default).
#'
#' @name renameClasses.SCM
#' @rdname renameClasses.SCM
#' @param class character; Old class names.
#' @return A character vector with new class names.
#' @author Jason Serviss
#'
NULL
#' @export
#' @importFrom dplyr case_when

renameClasses.SCM <- function(class) {
  case_when(
    class == "0" ~ "A375",
    class == "1" ~ "HCT116",
    class == "2" ~ "HOS",
    TRUE ~ "error"
  )
}