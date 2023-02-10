# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(data.table))
sh(library(ggseg))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

# THIS MUST BE RUN AFTER TABLE 2 HAS BEEN RUN

PATH.ABAGEN.EXPRESSION <- '../../tables/table2/abagen_expression_dkt.csv'
PATH.DKT.LABELS <- '../../tables/table2/dkt_labels.csv'
PATH.HITS <- '../../tables/table2/wightman_hits.csv'

for (path in c(PATH.ABAGEN.EXPRESSION, PATH.DKT.LABELS, PATH.HITS)) {
  if (! file.exists(path)) {
    msg = sprintf('Cannot find required input %s, make sure table2 is run in full.', path)
    stop(msg)
  }
}

# === read data ===========

expression <- fread(PATH.ABAGEN.EXPRESSION) # this takes a second
dkt.labels <- read.csv(PATH.DKT.LABELS) %>%
  filter(structure == 'cortex') %>%
  mutate(label = paste(tolower(hemisphere), 'h_', label, sep=''))
  # last bit converts to ggseg labels

# === load Wightman genes ==========

wightman.hits <- read.csv(PATH.HITS)

# === plot ============

gene.plot <- function(gene) {
  data = expression[, c('label', gene), with=F]
  plot.data <- left_join(dkt.labels, data, by=c('id'='label'))
  
  p <- ggplot(plot.data ) +
    geom_brain(atlas = dk,
               aes(fill=!!sym(gene))) +
    theme_void() + 
    scale_fill_colormap(colormap='magma') +
    theme(legend.position = 'none')
  
  return(p)
}


for (gene in unique(wightman.hits$Gene)) {
  print(gene)
  p <- gene.plot(gene)
  figname <- paste(gene, '.png', sep='')
  ggsave(figname, p, width=8, height=2)
}

# create final plot with big legend
p <- gene.plot(gene) +
  theme(legend.position='right',
        legend.key.size = unit(2, 'cm')) +
  guides(fill = guide_colorbar(ticks=F, label=F))
ggsave('BIG_LEGEND.png', p, width=8, height=8)

