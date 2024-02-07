
# imports
sh <- suppressPackageStartupMessages

sh(library(dplyr))
sh(library(stringr))
sh(library(this.path))
sh(library(R.matlab))

# set working directory to make relative paths work
setwd(this.dir())

# load necessary info for formatting
mat <- readMat('../../nmf/adni/results/mat/NumBases8.mat')
ptc.names <- read.csv('../../derivatives/adni/names_8ptc.csv')$Component
ptc.order <- read.csv('../../derivatives/adni/wscore_stage_order.csv')$Region
ptc.order <- str_replace(ptc.order, 'Cmp.', '')
region.order <- read.csv('../../derivatives/adni/nmf_regions_ggseg.csv')$label

# raw W matrix
W <- as.data.frame(mat$W)
colnames(W) <- ptc.names
W$Region <- region.order
W <- W %>%
  select(Region, all_of(ptc.order))
write.csv(W, 'ptcs.csv', row.names = F, quote = F, na = '')

# normalized to sum to 1
W <- as.data.frame(mat$Wnorm)
colnames(W) <- ptc.names
W$Region <- region.order
W <- W %>%
  select(Region, all_of(ptc.order))
write.csv(W, 'ptcs_norm.csv', row.names = F, quote = F, na = '')