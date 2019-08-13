renameClasses.MGA.abb <- function(class) {
  case_when(
    class == "0" ~ "C Stem 1",
    class == "1" ~ "C Goblet Plet1-",
    class == "2" ~ "C Colonocytes",
    class == "3" ~ "SI Stem",
    class == "4" ~ "SI Transient amplifying",
    class == "5" ~ "C Stem 3",
    class == "6" ~ "C Progenitor",
    class == "7" ~ "C Transient amplifying",
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
}
  
renameClasses.MGA <- function(class) {
  case_when(
    class == "0" ~ "Stem 1",
    class == "1" ~ "Goblet Plet1-",
    class == "2" ~ "Colonocytes",
    class == "3" ~ "Stem",
    class == "4" ~ "Transient amplifying",
    class == "5" ~ "Stem 3",
    class == "6" ~ "Progenitor",
    class == "7" ~ "Transient amplifying",
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

classOrder.MGA <- c(
  "SI.Paneth", "SI.Lgr5+", "SI.Lgr5+.Mki67+", "SI.TA.early", "SI.TA.late", "SI.Enterocytes", "SI.Goblet", 
  "C.Goblet.Mki67+", "C.Goblet.proximal", "C.Goblet.distal.Plet1", "C.Lgr5+.proximal.2",
  "C.Lgr5+.proximal.1", "C.Lgr5+.distal", "C.Lgr5+.Mki67+", "C.TA.distal", "C.Colonocytes", 
  "C.Goblet.distal.1",  "Enteroendocrine", "Tufft", "Blood"
)

renameClasses.SCM <- function(class) {
  case_when(
    class == "0" ~ "A375",
    class == "1" ~ "HCT116",
    class == "2" ~ "HOS",
    TRUE ~ "error"
  )
}