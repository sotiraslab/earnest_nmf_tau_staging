
# --- Imports ---------

library(dplyr)
library(this.path)

# --- Constants ---------

PTC.NAMES <- c('MedialTemporal',
               'LeftParietalTemporal',
               'RightParietalTemporal',
               'Precuneus',
               'Occipital',
               'LateralFrontal',
               'Sensorimotor',
               'Orbitofrontal')

# --- Data loading functions -------

load.ptc.matrix <- function(normalized = T) {
  
  # Overview
  # --------
  #
  # Load the spatial distribution of each PTC as a dataframe
  # This is the "W" matrix from NMF
  # These are saved in the "ptcs" folder of this directory
  
  # Parameters
  # ----------
  # normalized (logical): load W-matrix normalized so that individual PTCs
  #     sum to 1.
  
  target <- if(normalized) 'ptcs_norm.csv' else 'ptcs.csv'
  path <- file.path(this.dir(), 'ptcs', target)
  ptcs <- read.csv(path)
  return (ptcs)
}

load.wscore.parameters <- function(model, ptc) {
  
  # Overview
  # --------
  #
  # Load the parameters of the repeated W-score models for a
  # given PTC.  These are saved in the "wscores" directory.
  # Each CSV contains the intercept, coefficients (age & sex),
  # and residual for models estimated in the manuscript.
  
  # Parameters
  # ----------
  # model ('adni', 'oasis', or 'both'): Sets which W-score models to load.
  #     Use 'adni' for the models estimated from ADNI-CU, and 'oasis' for
  #     the models estimated from OASIS-CU. 'both' will concatenate the
  #     two of these.
  # ptc (character) : Which PTC to load the models for.  See `PTC.NAMES`
  #     in this file for acceptable choices.
  
  model <- tolower(model)
   
  if (model %in% c('adni', 'oasis')) {
    path <- file.path(this.dir(), 'wscores', model, sprintf('%s.csv', ptc))
    params <- read.csv(path)
    return (params)
  } else if (model == 'both') {
    adni <- load.wscore.parameters('adni', ptc)
    oasis <- load.wscore.parameters('oasis', ptc)
    both <- rbind(adni, oasis)
    return (both)
  } else  {
    stop('`model` must be "adni", "oasis", or "both".')
  }
}

tau.roi.names <- function() {
  # Show the names of tau ROIs in order used in the PTC matrix.
  df <- load.ptc.matrix()
  return (df$Region)
}

# --- Helper functions --------

assign.stages <- function(data, regions, stage.grouping, p='any', atypical=NA) {
  
  # Overview
  # --------
  #
  # Function for converting binary assessments of pathology into stage
  # assignments.  Applies typical staging logic of finding the highest stage
  # such that a subject is positive for pathology in that stage and all
  # preceding ones.  Uses a different labeling for people who meet criteria
  # for a higher stage but not all preceding ones.
  
  # For PTC staging specifcally, this function is not needed.
  # See `ptc.staging()`, below, which calls this.
  
  # Parameters
  # ----------
  # data (data.frame) : Dataset (subjects x features)
  # regions (character) : Columns to stage on.  Each entry should
  #     correspond to the name of a binary column in `data`.
  # stage.grouping (integer) : A vector of integers used to group regions
  #     into stages. This should be the same length as `regions` and contain
  #     a non-decreasing collection of integers starting from 1
  #     E.g., `c(1, 1, 2, 3)` groups the first two regions into one stage,
  #     and puts the 3rd and 4th into individual stages.
  # p (float in [0, 1], 'any', or 'all') : Sets how many sub-regions a person
  #     has to have pathology in to meet criteria for a stage.  If "all",
  #     they must be positive in all regions of the stage.  If "any", they
  #     must be positive in 1 region for the stage.  Otherwise, can be 
  #     a float between 0 and 1 to signify the proportion of regions the
  #     person should be positive for.  E.g., use `0.5` to indicate that
  #     the person must show positivity in more than half of the regions
  #     corresponding to a given stage.  This is more useful when 
  #     having stages that include many sub-regions.
  # atypical (numeric or character) : Symbol to use for people who are
  #     atypical/non-stageable.
  
  staged.data <- data.frame(id=1:nrow(data))
  unique.stages <- unique(stage.grouping)
  n <- length(unique.stages)
  for (i in sort(unique.stages)) {
    regions.current <- as.character(regions[stage.grouping == i])
    sub.data <- data[, regions.current, drop=F]
    ps <- rowSums(sub.data) / length(regions.current)
    if (p == 'any') {
      positive.for.stage <- as.numeric(ps > 0)
    } else if (p == 'all') {
      positive.for.stage <- as.numeric(ps == 1)
    } else {
      positive.for.stage <- as.numeric(ps >= p)
    }
    staged.data[[sprintf('positive.%s', i)]] <- positive.for.stage
  }
  
  staged.data <- dplyr::select(staged.data, -id)
  
  diffs <- apply(staged.data, 1, diff)
  if (n == 2) {
    increasing <- diffs <= 0
  } else {
    increasing <- apply(diffs <= 0, 2, all)
  }
  
  stage <- ifelse(increasing, rowSums(staged.data), atypical)
  return(stage)
} 

# --- User functions ----------



ptc.uptake <- function(data, tau.roi.columns = NULL) {
  
  # Overview
  # --------
  #
  # Calculate the tau uptake in each PTC for a new dataset.  This is done
  # by taking the matrix multiplication of the normalized PTC matrix (W)
  # and a new dataset of regional tau uptakes.
  
  # Parameters
  # ----------
  # data (data.frame): Dataset (subjects x features)
  # tau.roi.columns (character) : Indicate the columns of `data` which
  #      are the SUVR uptakes in each GM region.  This should be the 
  #      same 68 cortical GM regions of the FreeSurfer atlas, in the same
  #      order as used in the manuscript.  To see this order, call
  #      `tau.roi.names()`.  The convention used is the same as ggseg.
  
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

ptc.wscores <- function(data, model = 'adni', age.column = 'Age', sex.column = 'Sex',
                        cutoff = NULL, ptc.columns = NULL) {
  # Overview
  # --------
  #
  # Apply the models learned in the manuscript to convert
  # tau SUVRs in PTCs into W-scores.  The parameters of these models
  # are saved in the "wscores" folder.  These are reloaded and used
  # to recontruct the linear models used to model the normative age- and sex-
  # predicted tau uptake.
  
  # The new data MUST include a column for age and a column for sex, in order
  # to apply the models estimated in the manuscript.
  
  # Parameters
  # ----------
  # data (data.frame): Dataset (subjects x features) containing uptakes in each
  #     PTC.
  # model ('adni', 'oasis', or 'both'): Which W-score models to use. 
  #     see `load.wscore.parameters()` for more information.
  # age.column (character): Column in `data` which contains age.
  # sex.column (character): Column in `data` which contains sex, as a binary
  #     variable (1 = male, 0 = female).
  # cutoff (numeric): Binarize the W-scores at the given cutoff.  If NULL,
  #     scores are left as continuous values.
  # ptc.columns (character): Columns which contain tau uptake in each PTC.
  #     If `NULL`, uses `PTC.NAMES`.
  
  if (is.null(ptc.columns)) {
    ptc.columns <- PTC.NAMES
  }
  
  wscores <- matrix(data = NA, nrow = nrow(data), ncol = 8)
  for (i in seq_along(ptc.columns)) {
    ptc <- ptc.columns[i]
    w.parameters <- load.wscore.parameters(model = model, ptc = ptc)
    
    # read observed age/sex
    observed.uptake <- data[[ptc]]
    observed.predictors <- data[, c(age.column, sex.column)]
    if (! all(sort(unique(data[[sex.column]])) == c(0, 1))) {
      stop('`sex.column` is not a binary column (male=1, female=0)')
    }
    
    # collect into matrix
    observed.mat <- matrix(NA, nrow=nrow(observed.predictors), ncol=3)
    observed.mat[, 1] = 1. # intercept
    observed.mat[, 2] = observed.predictors[, 1]
    observed.mat[, 3] = observed.predictors[, 2]
    
    # model    
    w.parameter.matrix <- t(as.matrix(w.parameters[, c('Intercept', 'Age', 'GenderMale')]))
    w.residual <- w.parameters[['Residual']]
    estimated <- observed.mat %*% w.parameter.matrix
    residuals <- observed.uptake - estimated
    # d'oh for rowwise division   https://stackoverflow.com/a/20596490/13386979
    relative.residuals <- t(t(residuals) / w.residual)
    wscore <- rowMeans(relative.residuals)
    
    # save
    wscores[, i] <- wscore
  }
  if (! is.null(cutoff)) {
    wscores <- ifelse(wscores >= cutoff, 1, 0)
  }
  
  wscores <- as.data.frame(wscores)
  colnames(wscores) <- PTC.NAMES
  return (wscores)
}

ptc.staging <- function(positivity, ptc.columns = NULL, non.stageable = 'NS') {
  
  # Overview
  # --------
  #
  # Takes a dataframe of binary tau uptakes in each PTC and returns the assigned
  # PTC stages.
  
  # Parameters
  # ----------
  # positivity (data.frame): Dataset (subjects x features) which contains
  #     columns with the binary assessment of tau pathology in each PTC.
  # ptc.columns (chatacter): If not using the default names of PTCs (`PTC.NAMES`),
  #     lets the user specify the 8 PTC columns.
  # non.stageable (character or numeric): Symbol to use to signify non-stageable
  #     individuals.
  
  if (is.null(ptc.columns)) {
    ptc.columns <- PTC.NAMES
  }
  
  stages <- assign.stages(data = positivity,
                          regions = ptc.columns,
                          stage.grouping = c(1, 2, 2, 2, 3, 3, 4, 4),
                          p = 'any', 
                          atypical = non.stageable)
  return (stages)
}

ptc.staging.pipeline <- function(data,
                                 tau.roi.columns = NULL,
                                 model = 'adni',
                                 age.column = 'Age',
                                 sex.column = 'Sex',
                                 cutoff = 2.5,
                                 non.stageable = 'NS') {
  
  # Overview
  # --------
  #
  # Stitches all the PTC staging steps into 1 pipeline.  See the following
  # functions for more detail: ptc.uptake(), ptc.wscores(), and ptc.staging().
  
  # Parameters
  # ----------
  # data (data.frame): Dataset (subjects x features)
  # tau.roi.columns (character) : Indicate the columns of `data` which
  #      are the SUVR uptakes in each GM region.  This should be the 
  #      same 68 cortical GM regions of the FreeSurfer atlas, in the same
  #      order as used in the manuscript.  To see this order, call
  #      `tau.roi.names()`.  The convention used is the same as ggseg.
  # model ('adni', 'oasis', or 'both'): Which W-score models to use. 
  #     see `load.wscore.parameters()` for more information.
  # age.column (character): Column in `data` which contains age.
  # sex.column (character): Column in `data` which contains sex, as a binary
  #     variable (1 = male, 0 = female).
  # cutoff (numeric) : Cutoff to use when binarizing W-scores.
  # non.stageable (character or numeric): Symbol to use to signify non-stageable
  #     individuals.
  #
  # Returns
  # -------
  # output (list): a list containing the PTC uptakes, W-scores, binarized
  #     W-scores, and stages for the intput dataset.
  
  uptakes <- ptc.uptake(data, tau.roi.columns = tau.roi.columns)
  uptakes.big <- cbind(uptakes, data[, c(age.column, sex.column)])
  wscores <- ptc.wscores(uptakes.big,
                         model = model,
                         age.column = age.column,
                         sex.column = sex.column,
                         cutoff = NULL)
  positivity <- as.data.frame(ifelse(wscores >= cutoff, 1, 0))
  stages <- ptc.staging(positivity, non.stageable = non.stageable)
  
  output <- list(
    'uptakes' = uptakes,
    'wscores' = wscores,
    'positivity' = positivity,
    'stages' = stages
  )
  
  return (output)
}
