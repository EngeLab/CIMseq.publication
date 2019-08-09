library(CIMseq.testing)
library(tidyverse)
library(ggdendro)
library(patchwork)

###singlets
singlets.1 <- syntheticSinglets(5000, 20, 1, seed = 74283)
singlets.2 <- syntheticSinglets(5000, 30, 1, seed = 89323)
singlets.3 <- syntheticSinglets(5000, 10, 1, seed = 23934)
singlets.4 <- syntheticSinglets(5000, 25, 1, seed = 87923)

colnames(singlets.1) <- paste0(colnames(singlets.1), ".", 1)
colnames(singlets.2) <- paste0(colnames(singlets.2), ".", 2)
colnames(singlets.3) <- paste0(colnames(singlets.3), ".", 3)
colnames(singlets.4) <- paste0(colnames(singlets.4), ".", 4)

singlets <- do.call(cbind, list(singlets.1, singlets.2, singlets.3, singlets.4))
rownames(singlets) <- paste0("gene", 1:nrow(singlets))

fc <- foldChangePerClass(
  singlets, 
  tibble(
    sample = colnames(singlets), 
    class = str_replace(colnames(singlets), ".*\\.(.*)", "\\1")
  )
)

features <- apply(fc, 1, function(r) length(which(r > 4)) == 1)
features <- which(features)
names(features) <- apply(fc[features, ], 1, function(r) colnames(fc)[which(r > 4)])

singlets.cpm <- CIMseq:::.norm.counts(singlets)
singlets.log <- CIMseq:::.norm.log.counts(singlets)

###multiplets
#desired combos: 2-1, 2-4, random
nMultipletsPerClass <- 30
set.seed(5678)
fractions.1 <- map(1:nMultipletsPerClass, function(i) {
  f <- c(runif(2, 0.25, 0.8), runif(2, 0, 0.01))
  f / sum(f)
})
fractions.2 <- map(1:nMultipletsPerClass, function(i) {
  f <- c(runif(2, 0.25, 0.8), runif(2, 0, 0.01))
  f <- f[c(3, 1, 4, 2)]
  f / sum(f)
})
fractions.3 <- map(1:nMultipletsPerClass, function(i) {
  f <- c(runif(2, 0.25, 0.8), runif(6, 0, 0.01))
  f <- sample(f, 4, FALSE)
  f / sum(f)
})

# fractions.3 <- list(
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0.25, 0.8), runif(1, 0, 0.01)),
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0.25, 0.8), runif(1, 0, 0.01)),
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0.25, 0.8), runif(1, 0, 0.01)),
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0, 0.01), runif(1, 0.25, 0.8)),
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0, 0.01), runif(1, 0.25, 0.8)),
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0, 0.01), runif(1, 0.25, 0.8)),
#   c(runif(1, 0.25, 0.8), runif(1, 0, 0.01), runif(1, 0, 0.01), runif(1, 0.25, 0.8)),
#   c(runif(1, 0.25, 0.8), runif(3, 0, 0.01)),
#   c(runif(1, 0.25, 0.8), runif(3, 0, 0.01)),
#   c(runif(3, 0, 0.01), runif(1, 0.25, 0.8))
# )
# fractions.3 <- map(fractions.3, function(f) f / sum(f))
fractions <- c(fractions.1, fractions.2, fractions.3)

baseSeed <- 2389
multiplets <- map2(1:(nMultipletsPerClass * 3), fractions, function(i, f) {
  syntheticMultipletsFromCounts(
    singlets.cpm,
    str_replace(colnames(singlets), ".*\\.(.*)", "\\1"),
    f,
    baseSeed + i
  )
}) %>%
  reduce(cbind)

colnames(multiplets) <- c(
  paste0(rep("2-1", nMultipletsPerClass), ".", 1:nMultipletsPerClass),
  paste0(rep("2-4", nMultipletsPerClass), ".", 1:nMultipletsPerClass),
  paste0(rep("r", nMultipletsPerClass), ".", 1:nMultipletsPerClass)
)

####PLOT
#heat function and data
classHeatmap <- function (data, counts.log, classes, class.order) 
{
  geneOrd <- function(c, data, counts.log) {
    IDs <- filter(data, class == c)$gene
    if (length(IDs) == 1) {
      IDs
    }
    else {
      mat <- counts.log[rownames(counts.log) %in% IDs, ]
      cr <- cor(t(mat), method = "pearson")
      my.dist <- as.dist(1 - cr)
      c <- hclust(my.dist)
      c[[4]][c[[3]]]
    }
  }
  
  order <- map(unique(feature.df$class), ~geneOrd(.x, feature.df, singlets.log)) %>% 
    setNames(unique(feature.df$class)) %>% 
    namedListToTibble(.) %>% 
    rename(class = names, gene = variables) %>% 
    mutate(order = 1:n())
  
  p <- counts.log[rownames(counts.log) %in% data$gene, colnames(counts.log) %in% classes$sample] %>% 
    apply(., 1, function(x) (x - min(x))/(max(x) - min(x))) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene") %>% 
    as_tibble() %>% 
    gather(sample, value, -gene) %>% 
    left_join(classes, by = "sample") %>% 
    rename(sampleClass = class) %>% 
    arrange(sampleClass) %>% 
    mutate(sampleClass = parse_factor(sampleClass, levels = class.order)) %>%
    mutate(plotSample = factor(sample, levels = unique(sample))) %>% 
    left_join(order, by = "gene") %>% 
    rename(geneClass = class) %>% 
    arrange(order) %>% 
    mutate(plotGene = factor(gene, levels = unique(gene))) %>%
    ggplot(aes(plotSample, plotGene)) + 
    geom_tile(aes(fill = value)) + 
    facet_grid(. ~ sampleClass, scales = "free", space = "free_x") + 
    scale_fill_viridis() + 
    theme_few() + 
    theme(
      legend.position = "bottom", 
      axis.text.y = element_blank(), axis.text.x = element_blank(), 
      axis.ticks = element_blank(), strip.text.x = element_text(angle = 90), 
      strip.text.y = element_text(angle = 0)
    ) + 
    labs(x = "Samples", y = "Genes") + 
    guides(fill = guide_colourbar(
      title = "z-score", title.position = "top", 
      title.hjust = 0.5, barwidth = 10
    ))
  p
  return(p)
}

feature.df <- data.frame(
  gene = rownames(singlets)[features],
  class = names(features),
  stringsAsFactors = FALSE
)

#stacked bar
bar <- fractions %>%
  reduce(rbind) %>%
  magrittr::set_colnames(sort(unique(names(features)))) %>%
  magrittr::set_rownames(colnames(multiplets)) %>%
  matrix_to_tibble("sample") %>%
  gather(class, frac, -sample) %>%
  ggplot() +
  geom_bar(aes(sample, frac, fill = class), stat = "identity") +
  ggthemes::scale_fill_ptol() +
  ggthemes::theme_few() +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank()
  ) +
  labs(y = "Fraction") +
  guides(fill = guide_legend(title = NULL))


#dendrogram
h.clust <- hclust(dist(t(singlets)))
cOrder <- unique(str_replace(h.clust[[4]][h.clust[[3]]], ".*\\.(.*)", "\\1"))
dendro <- ggdendrogram(h.clust) + theme_dendro()

#heat
heat <- classHeatmap(
  feature.df,
  cbind(singlets.cpm, multiplets),
  tibble(
    sample = c(colnames(singlets), colnames(multiplets)),
    class = c(
      str_replace(colnames(singlets), ".*\\.(.*)", "\\1"),
      str_replace(colnames(multiplets), "(.*)\\..*", "\\1")
    )
  ),
  c(cOrder, "2-1", "2-4", "r")
)

#all
dendro + bar - heat + plot_layout(ncol = 1)

#circos
.closestCircle <- function(from, to, class, max) {
  o <- class == to
  to[o] <- from[o]
  from[o] <- class[o]
  if(any(is.na(from))) {
    cat(paste("NA pos: ", from))
  }
  mid <- as.integer(max/2)
  adj <- mid-from
  from <- (from+adj) %% max
  to <- (to+adj) %% max
  pos <- to-from
  pos[pos < 0] <- -mid-pos[pos < 0]
  pos[pos > 0] <- mid-pos[pos > 0]
  return(pos)
}

#2 = green, 1 = blue, 3 = yellow, 4 = pink
colours <- tibble(
  class = cOrder, 
  colour = c("#117733", "#DDCC77", "#4477AA", "#CC6677"),
  nr = 1:length(cOrder),
  combined = paste0("(", nr, ") ", class)
)

connections <- map(fractions, function(f) {
  ct <- as.character(1:4)[f > 0.25]
  if(length(ct) == 1) {
    NA
  } else {
    apply(combn(ct, 2), 2, function(c) paste(sort(c), collapse = "_"))
  }
})
names(connections) <- paste0("s", 1:length(connections))
connections <- unlist(connections[!is.na(connections)])

data <- tibble(
  sample = names(connections), 
  connection = connections, 
  connectionID = 1:length(connections)
) %>% 
  mutate(
    from = map_chr(connection, ~str_split(.x, "_")[[1]][1]),
    to = map_chr(connection, ~str_split(.x, "_")[[1]][2])
  ) %>% 
  gather(direction, class, -sample, -connection, -connectionID) %>% 
  select(-direction) %>% 
  mutate(
    from = map_chr(connection, ~str_split(.x, "_")[[1]][1]),
    to = map_chr(connection, ~str_split(.x, "_")[[1]][2])
  ) %>%
  group_by(class) %>%
  arrange(.closestCircle(
    match(from, cOrder), match(to, cOrder), 
    match(class, cOrder), 12), 
    .by_group = TRUE
  ) %>%
  mutate(position = 1:n()) %>%
  ungroup() %>%
  mutate(significant = if_else(connection %in% c("1_2", "2_4"), TRUE, FALSE)) %>%
  mutate(p.col = if_else(significant, "#4B0B5A", "grey95")) %>%
  arrange(significant)

gap.degree <- 200 / length(cOrder)

posdat <- data %>% 
  group_by(class) %>% 
  summarize(start = min(position), end = max(position)) %>% 
  arrange(match(class, colours$class)) %>%
  as.data.frame() %>%
  column_to_rownames("class") %>%
  as.matrix()

circos.par(gap.degree = gap.degree, cell.padding = c(0, 0))
circos.initialize(factors = as.character(cOrder), xlim = posdat)
circos.trackPlotRegion(
  ylim = c(0, 1), bg.col = pull(colours, colour),
  bg.border = NA, track.height = 0.1
)

# for(i in 1:nrow(colours)) {
#   circos.text(
#     x = mean(posdat[as.character(pull(colours, class)[i]), ]), y = 0.5,
#     labels = as.character(pull(colours, class)[i]),
#     sector.index = pull(colours, class)[i], track.index = 1, col = "white",
#     facing = "downward", cex = label.cex
#   )
# }

#add links
for(i in 1:length(unique(data$connectionID))) {
  conn <- filter(data, connectionID == unique(data$connectionID)[i])
  if(nrow(conn) != 2) stop("error")
  circos.link(
    pull(conn, class)[1], pull(conn, position)[1],
    pull(conn, class)[2], pull(conn, position)[2],
    col = unique(pull(conn, p.col)), lwd = 2
  )
}

circos.clear()

#legend 
pal <- colorRampPalette(c("grey95", viridis::viridis(1)))(120)[30:120]
p <- ggplot() + 
  geom_point(aes(1:length(pal), 1:length(pal), colour = 1:length(pal))) + 
  scale_colour_gradientn(colours = c(pal[1], pal[length(pal)])) +
  guides(colour = guide_colorbar(
    title = "Obs. / Exp.", title.position = "top", title.hjust = 0.5
  )) +
  theme(
    legend.position = "top",
    legend.margin = margin(t = 0, b = 0, unit='cm'),
    legend.key = element_blank(),
    legend.box.background = element_blank()
  )
l <- g_legend(p)
draw_legend(l)

####SAVE EVERYTHING
