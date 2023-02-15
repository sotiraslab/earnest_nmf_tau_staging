# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(tidyverse))
sh(library(this.path))

# === Set WD ===========

setwd(this.dir())

# === Required Files ===========

PATH.ADNI <- '../../derivatives/adni/main_data_with_8ptc.csv'
PATH.OASIS <- '../../derivatives/oasis3/main_data_with_8ptc.csv'

# === Load Data ========

adni <- read.csv(PATH.ADNI)
oasis <- read.csv(PATH.OASIS)

# === Select components for analysis ==========

cols <- colnames(adni)
components <- cols[grepl('Cmp', cols)]

# === ADNI Normalization ========

control.vals <- adni %>%
  filter(Group == 'ControlBaseline') %>%
  select(all_of(components))

normalizer <- data.frame(Component=colnames(control.vals),
                         Mean=colMeans(control.vals),
                         SD=apply(control.vals, 2, sd))
rownames(normalizer) <- NULL

adni.normalized <- adni
for (i in 1:nrow(normalizer)) {
  row <- normalizer[i, ]
  cmp <- row$Component
  mu <- row$Mean
  sigma <- row$SD

  adni.normalized[[cmp]] <- (adni.normalized[[cmp]] - mu) / sigma
}

path.out <- '../../derivatives/adni/data_loadings_zscore.csv'
write.csv(adni.normalized, path.out, row.names=F, na='', quote=F)


# === same for OASIS ==========

control.vals <- oasis %>%
  filter(Group == 'ControlSet') %>%
  select(all_of(components))

normalizer <- data.frame(Component=colnames(control.vals),
                         Mean=colMeans(control.vals),
                         SD=apply(control.vals, 2, sd))
rownames(normalizer) <- NULL

oasis.normalized <- oasis
for (i in 1:nrow(normalizer)) {
  row <- normalizer[i, ]
  cmp <- row$Component
  mu <- row$Mean
  sigma <- row$SD
  
  oasis.normalized[[cmp]] <- (oasis.normalized[[cmp]] - mu) / sigma
}

path.out <- '../../derivatives/oasis3/data_loadings_zscore.csv'
write.csv(oasis.normalized, path.out, row.names=F, na='', quote=F)
