#' plotPDF
#'
#' Plotting theme
#'
#' @name plotPDF
#' @rdname plotPDF
#' @param plot A ggplot
#' @param strip.text.size Numeric size.
#' @return A ggplot.
#' @author Jason T. Serviss

#' @export
#' @importFrom ggthemes scale_fill_ptol scale_colour_ptol

plotPDF <- function(plot, strip.text.size = 3.25) {
  p <- plot +
    theme_bw() +
    scale_fill_ptol() +
    scale_colour_ptol() +
    theme(
        plot.title = element_blank(),
        plot.margin = unit(rep(0.1, 4), "lines"),
        legend.position = "top",
        legend.background = element_blank(),
        legend.box.spacing = unit(5, "pt"), #this controls the spacing between strip.text.x and the legend
        legend.margin = margin(rep(0, 4), unit = "pt"),
        legend.key.size = unit(1/5, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.length = unit(1/15, "cm"),
        strip.text.x = element_text(
          margin = margin(0.01, 0, 5, 0, "pt")
        ),
        strip.text = element_text(size = strip.text.size),
        panel.border = element_rect(fill = NA, size = 0.15)
    )
    return(p)
}

#' plotRmarkdown
#'
#' Plotting theme
#'
#' @name plotRmarkdown
#' @rdname plotRmarkdown
#' @param plot A ggplot
#' @return A ggplot.
#' @author Jason T. Serviss

#' @export
#' @importFrom ggthemes scale_fill_ptol scale_colour_ptol

plotRmarkdown <- function(plot) {
  p <-  plot +
    theme_bw() +
    scale_colour_ptol() +
    scale_fill_ptol() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 15),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 20
      ),
      plot.caption = element_text(
        hjust = 0,
        margin = margin(t = 15),
        family = "Times New Roman",
        size = 14
      ),
      axis.title = element_text(size = 17),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.text = element_text(size = 12)
    )
  
  return(p)
}
