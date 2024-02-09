
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
  target <- if(normalized) 'ptcs_norm.csv' else 'ptcs.csv'
  path <- file.path(this.dir(), 'ptcs', target)
  ptcs <- read.csv(path)
  return (ptcs)
}

load.wscore.parameters <- function(model, ptc) {
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

# --- Helper functions --------

assign.stages <- function(data, regions, stage.grouping, p='any', atypical=NA) {
  
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
