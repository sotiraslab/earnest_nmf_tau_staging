
# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(anticlust))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Load ==========

df <- read.csv("../../derivatives/oasis3/main_data.csv") %>%
  filter(Group == 'TrainingSet')

# === Load NMF Creation function ==========

source("../../scripts/create_nmf_input.R")

# === Create NMF input and save outputs ==========

features = colnames(df)[grepl('CTX', colnames(df))]
outdata = '../../derivatives/oasis3/nmf_matrix.csv'
outfeatures = '../../derivatives/oasis3/nmf_regions.csv'
outids = '../../derivatives/oasis3/nmf_ids.csv'

X <- create.nmf.input(df, features=features, id.cols = c("Subject", "Session"),
                      outdata = outdata, outfeatures = outfeatures,
                      outids = outids)

# ==== Repeated reproducibility splits =========

set.seed(42)

# number of repeats
N <- 50

splits <- matrix(data=NA, nrow=nrow(df), ncol=N)
continuous.vars <- df[, c("Age", "TotalCtxTauMean")]
categorical.vars <- df[, c("GENDER", "CDR")]

for (i in 1:N) {
  v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
  splits[, i] <- v
}

splits <- as.data.frame(splits)
outpath <- '../../derivatives/oasis3/reproducibility_split_indices.csv'
write.table(splits, outpath, row.names = F, col.names = F, sep=',')
