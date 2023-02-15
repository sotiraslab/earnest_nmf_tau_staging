# === Import =======

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggseg))
sh(library(R.matlab))
sh(library(this.path))
sh(library(tidyverse))

# === Set working directory =======

setwd(this.dir())

# === Required Files =======

PATH.PTC8.MAT <- '../../nmf/adni/results/mat/NumBases8.mat'
PATH.NAMES.PTC8 <- '../../derivatives/adni/names_8ptc.csv'
PATH.REGIONS <- '../../derivatives/adni/nmf_regions_ggseg.csv'

# === Define WTA function ==========

winner.take.all.assignment <- function(W, norm=F) {
  
  if (norm) W <- W / apply(W, 2, sum)
  assignment <- apply(W, 1, which.max)
  return (assignment)
  
}

# === Read in components ========

mat <- readMat(PATH.PTC8.MAT)
W <- mat$Wnorm

# === Assignment =========

assignment <- winner.take.all.assignment(W)

# === Create output dataframe =========

cmp.cols <- read.csv(PATH.NAMES.PTC8)$Component
cmp.cols.long <- paste('Cmp.', cmp.cols, sep='')
mapper <- data.frame(value=1:8, name=cmp.cols.long)


# read in specific labels to make ggseg compatible
regions <- read.csv(PATH.REGIONS)

assignment.df <- regions
assignment.df$value <- assignment
assignment.df <- left_join(assignment.df, mapper, by='value')
assignment.df$value <- factor(assignment.df$value)

path.out <- '../../derivatives/adni/ptc8_winner_take_all.csv'
write.csv(assignment.df, path.out)

# === Plot ==========


colors = colormap('jet', nshades = 8)
color.mapper <- cmp.cols
names(color.mapper) <- as.character(1:8)

ggplot(assignment.df) +
  geom_brain(atlas = dk,
             color='black',
             aes(fill=value),
             position = position_brain(hemi ~ side )) +
  theme_void() +
  scale_fill_manual(values=colors, labels=cmp.cols)

ggsave('assignment_figure.png', width=8, height=8)

