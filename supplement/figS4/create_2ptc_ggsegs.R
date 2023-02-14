
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
PATH.NMF.ADNI <- '../../nmf/adni/results/mat/NumBases2.mat'
PATH.NMF.OASIS <- '../../nmf/oasis3/results/mat/NumBases2.mat'
PATH.MATCH <- '../figS1/adni_v_oasis_compare/matching/Match2.mat'

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
               color='black',
               position = position_brain(hemi ~ side)) +
    theme_void() +
    theme(legend.position = 'none',
          plot.margin=unit(c(0, 0, 0, 0),"cm")) + 
    scale_fill_colormap(colormap='viridis')
    
  
  return (p)
}

# === Read region labels for ggseg ==========

regions.df <- read.csv(PATH.REGIONS)
REGION.LABELS <- regions.df$label

# === Read matching ========

matches <- readMat(PATH.MATCH)
match.idx <- matches$idx.hug1

# === plot ADNI ==========

for (i in 1:2) {
  p <- plot_ptc_ggseg(results.mat = PATH.NMF.ADNI,
                      N = i,
                      region.labels = REGION.LABELS)
  savename <- sprintf('adni_ptc_%s.png', i)
  ggsave(savename, p, width=7, height=7, units='in')
}

# === plot OASIS ==========

for (i in 1:2) {
  idx.oasis <- match.idx[i]
  p <- plot_ptc_ggseg(results.mat = PATH.NMF.OASIS,
                 N = idx.oasis,
                 region.labels = REGION.LABELS)
  savename <- sprintf('oasis_ptc_%s.png', i)
  ggsave(savename, p, width=7, height=7, units='in')
}
