
# --- Imports ---------

library(dplyr)
library(this.path)

# --- Data loading functions -------

load.ptc.matrix <- function(normalized = T) {
  target <- if(normalized) 'ptcs_norm.csv' else 'ptcs.csv'
  path <- file.path(this.dir(), 'ptcs', target)
  ptcs <- read.csv(path)
  return (ptcs)
}

# --- Helper functions --------

# --- User functions ----------

ptc.uptake <- function(data, tau.roi.columns = NULL) {
  ptcs <- load.ptc.matrix(normalized = T)
  if (is.null(tau.roi.columns)) {
    tau.roi.columns <- ptcs$Region
  }
  X <- as.matrix(data[, tau.roi.columns]) #  n x 68
  W <- as.matrix(select(ptcs, -Region))   # 68 x  8
  uptakes <- X %*% W                      #  n x  8
  uptakes.df <- as.data.frame(uptakes)
  colnames(uptakes.df) <- colnames(W)
  return (uptakes.df)
}

ptc.wscores <- function(data, model = 'ADNI', age.column = 'Age', sex.column = 'Sex',
                        cutoff = NULL, ptc.columns = NULL) {
  
}

ptc.staging.pipeline <- function(tau.rois) {
  
}

df <- read.csv('../../derivatives/adni/main_data.csv')
rois <- colnames(df)[grepl('CTX_.*_SUVR', colnames(df), perl = T)]
uptakes <- ptc.uptake(df, tau.roi.columns = rois)