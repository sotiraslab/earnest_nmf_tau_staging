# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(R.matlab))
sh(library(this.path))
sh(library(tidyverse))

# === Set wd ===========

setwd(this.dir())


# === Read full dataset ==========

df.path <- '../../derivatives/adni/main_data.csv'
df <- read.csv(df.path)

# === read NMF stuff ==========

mat <- readMat('../../nmf/adni/results/mat/NumBases8.mat')
regions <- read.csv('../../derivatives/adni/nmf_regions.csv')$Feature

# === Project ===========

Wnorm <- mat$Wnorm
Xnew <- as.matrix(df[, regions])
projected <- Xnew %*% Wnorm

# === Reassign to dataframe ==========

comp.names <-  c(
  "Cmp.Orbitofrontal",
  "Cmp.LateralFrontal",
  "Cmp.MedialTemporal",
  "Cmp.Precuneus",
  "Cmp.LeftParietalTemporal",
  "Cmp.Sensorimotor",
  "Cmp.Occipital", 
  "Cmp.RightParietalTemporal"
)
projected <- as.data.frame(projected)
colnames(projected) <- comp.names
df <- cbind(df, projected)

# === Save outputs ==========

outpath <- '../../derivatives/adni/main_data_with_8ptc.csv'
write.csv(df, outpath, na='', quote=F, row.names = F)

drop.rois <- df %>%
  select(-contains('CTX'))

outpath <- '../../derivatives/adni/main_data_with_8ptc_no_rois.csv'
write.csv(drop.rois, outpath, na='', quote=F, row.names = F)

# component names, as ordered by NMF output
temp <- data.frame(Component=gsub('Cmp.', '', comp.names))
write.csv(temp, '../../derivatives/adni/names_8ptc.csv', row.names = F)