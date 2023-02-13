
# === Import =======

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggplot2))
sh(library(ggseg))
sh(library(R.matlab))
sh(library(stringr))
sh(library(this.path))

# === Set WD ==========

setwd(this.dir())

# === Files needed =========

PATH.REGIONS <- '../../derivatives/adni/nmf_regions_ggseg.csv'
PATH.PTC.NAMES <- '../../derivatives/adni/names_8ptc.csv'
PATH.NMF.MAT <- '../../nmf/adni/results/mat/NumBases8.mat'

# === Define plotting function ========

plot_ptc_ggseg <- function(results.mat, N, region.labels) {
  
  # === Read Components =========
  data <- readMat(results.mat)
  W <- data$W
  dims <- dim(W)
  FEATURES<- dims[1]
  RANK <- dims[2]
  
  
  # === Create plot data ==========
  PLOT.DATA <- data.frame(label = region.labels,
                          value = W[, N])
  
  # === Plot ==========
  p <- ggplot(PLOT.DATA) +
    geom_brain(atlas = dk,
               aes(fill=value),
               color='black') +
    theme_void() +
    theme(legend.position = 'none',
          plot.margin=unit(c(0, 0, 0, 0),"cm")) + 
    scale_fill_colormap(colormap='viridis')
  
  return (p)
}

# === Read region labels for ggseg ==========

regions.df <- read.csv(PATH.REGIONS)
REGION.LABELS <- regions.df$label

# === Read component names ========

cmp.names <- read.csv(PATH.PTC.NAMES)$Component

# === Run ==========

for (i in 1:8) {
  p <- plot_ptc_ggseg(results.mat = PATH.NMF.MAT,
                 N = i,
                 region.labels = REGION.LABELS)
  savename <- sprintf('%s.png', cmp.names[i])
  ggsave(savename, p, width=1800, height=1800/8, units='px')
}
