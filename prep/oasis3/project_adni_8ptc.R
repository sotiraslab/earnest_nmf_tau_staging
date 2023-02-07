# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(R.matlab))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Read ADNI components and labels ==========

path.components <- '../../nmf/adni/results/mat/NumBases8.mat'
path.regions <- '../../derivatives/adni/nmf_regions.csv'

mat <- readMat(path.components)
components <- mat$Wnorm # already normalized
regions <- read.csv(path.regions)$Feature

# === Read ADRC data =========

path.adrc <- '../../derivatives/oasis3/main_data.csv'
adrc <- read.csv(path.adrc)

# === Pull ADRC regional data ==========

# The following should grab the ROIs in the same order as is
# encoded by the ADNI components

adrc.suvrs <- adrc[, regions]

# === Project! =========

adrc.project <- as.matrix(adrc.suvrs) %*% as.matrix(components)

# === Add ADNI component names ==========

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
colnames(adrc.project) <- comp.names

# === Merge into ADRC DF ============

adrc.big <- cbind(adrc, as.data.frame(adrc.project))
adrc.edit <- dplyr::select(adrc.big, -all_of(regions))

# === Save ===========

outpath <- '../../derivatives/oasis3/main_data_with_8ptc_no_rois.csv'
write.csv(adrc.edit, outpath,
          na='', quote=F, row.names=F)

outpath <- '../../derivatives/oasis3/main_data_with_8ptc.csv'
write.csv(adrc.big, outpath, na='', quote=F, row.names=F)