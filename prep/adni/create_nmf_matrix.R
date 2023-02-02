
# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(anticlust))
sh(library(dplyr))
sh(library(lubridate))
sh(library(tidyr))
sh(library(this.path))

# === Set WD ==========

setwd(this.dir())

# === Load data for NMF ==========

PATH <- '../../derivatives/adni/main_data.csv'
df <- read.csv(PATH) %>%
  mutate(DateTau = as_datetime(ymd(DateTau))) %>%
  filter(Group == 'TrainingBaseline') 

# === Load helper func ==========

source("../../scripts/create_nmf_input.R")

# === prep NMF ==========

features = colnames(df)[grepl('CTX_.*_SUVR', colnames(df), perl=T)]
outdata = '../../derivatives/adni/nmf_matrix.csv'
outfeatures = '../../derivatives/adni/nmf_regions.csv'
outids = '../../derivatives/adni/nmf_ids.csv'

X <- create.nmf.input(df, features=features, id.cols = c("RID", "DateTau"),
                      outdata = outdata, outfeatures = outfeatures,
                      outids = outids)

# === Create reproducibility splits - repeated =============

set.seed(42)

# number of repeats
N <- 50

splits <- matrix(data=NA, nrow=nrow(df), ncol=N)
continuous.vars <- df[, c("Age", "CorticalTauAverage")]
categorical.vars <- df[, c("Gender", "CDRBinned")]

for (i in 1:N) {
  v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
  splits[, i] <- v
}

splits <- as.data.frame(splits)
outpath <- '../../derivatives/adni/reproducibility_split_indices.csv'
write.table(splits, outpath, row.names = F, col.names = F, sep=',')