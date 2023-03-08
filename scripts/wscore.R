
library(anticlust)

repeated.wscore.train <- function(control.data, y, covariates,
                                  match.continuous, match.categorical,
                                  portion.train = 0.8, repeats = 200,
                                  seed = T, verbose = T) {
  
  # Overview
  # --------
  
  # Train a repeated sampling W-score model. This is the basic method:
  
  #   1. Select a normative dataset relative to your population of interest.
  #   2. Split this dataset into "train" and "test" portions.
  #   3. In the training set, learn a linear model to predict a variable of
  #      interest from a set of covariates.
  #   4. In the test set, predict the outcome of interest and record the
  #      standard deviation of the residuals (normative residual).
  #   5. Repeat steps 2-4 a large number of times, saving each model and 
  #      normative residual.
  #   6. For predicting new data, predict on it with all models.  For each,
  #      get an individual w-score which is the model residual divided by
  #      the corresponding normative residual.  To get a final W-score,
  #      average all the individual w-scores from the repeats.
  #      (see repeated.wscore.predict)
  
  # This is based on an approach described in Lee et al. (Neuron, 2022).
  # The major difference is the splitting of the control set and the 
  # repeated resampling. The additional splitting of the control set is applied
  # to provide a more conservative estimate of the residuals.  Because of variance
  # in this splitting, repeated resampling and averaging is applied.

  
  # Parameters
  # ---------
  # control.data: (data.frame) Table containing wide-form control data
  # y (character): Character vector containing one or more columns to
  #     in the control.data to model W-scoring for
  # covariates (character): columns of the control.data to use as normative
  #     predictors in the linear models
  # match.continuous / match.categorical (character): Both character vectors
  #     specifying variables to use for matching control samples when splitting,
  #     supplied to anticlust.  Right now this is mandatory.
  # portion.train (float): proportion of training/testing when splitting controls
  # repeats (int): number of times to resample
  # seed (logical): When true, set the random seed each resampling.  Sets to the
  #     index of each sample.
  # verbose (logical): add some print statements.
  
  # returned at top level
  big.output <- list()
  return.splits <- matrix(NA, nrow=nrow(control.data), ncol=repeats)
  
  # first loop, creating data splits
  ntrain <- round(portion.train * nrow(control.data))
  ntest <- nrow(control.data) - ntrain
  
  control.train.holder <- list()
  control.test.holder <- list()
  
  for (i in 1:repeats) {
    if (verbose) cat("Data splitting repetition ", i, '/', repeats, "\n")
    
    continuous <- control.data[,  match.continuous, drop=F]
    categorical <- control.data[, match.categorical, drop=F]
    
    if (seed) set.seed(i)
    splits <- anticlustering(continuous, K=c(ntrain, ntest), categories = categorical)
    control.train <- control.data[splits == 1, ]
    control.test <- control.data[splits == 2,]
    control.train.holder[[i]] <- control.train
    control.test.holder[[i]] <- control.test
    return.splits[, i] <- splits
  }
  
  # second loop, create models
  for (j in 1:length(y)) {
    col <- y[j]
    if (verbose) cat("W scoring for y =", col, "\n")
    
    holder.models <- list()
    holder.residuals <- rep(NA, repeats)
    
    for (i in 1:repeats) {
      control.train <- control.train.holder[[i]]
      control.test <- control.test.holder[[i]]
      
      fml <- paste(col, '~', paste(covariates, collapse = ' + '), collapse=' ')
      model <- lm(fml, data=control.train)
      holder.models[[i]] <- model
      
      predicts.control <- predict(model, control.test)
      residuals.control <- control.test[[col]] - predicts.control
      residual.control.sd <- sd(residuals.control)
      holder.residuals[[i]] <- residual.control.sd
      
    }
    collect = list(
      models = holder.models,
      residuals = holder.residuals
    )
    big.output[[col]] <- collect
  }
  
  return (big.output)
  
}

repeated.wscore.predict <- function(output, new.data, predict.vars = NULL) {
  
  # Overview
  # --------
  
  # Apply a trained W-score model on new data.  See repeated.wscore.train
  # for modeling description and model training.
  
  # Parameters
  # ---------
  # output (list): output of `repeated.wscore.train`
  # new.data (data.frame): table to predict new W-scores for
  # predict.vars (character): specific columns to predict.  When not specified,
  #     tries to predict for entries in `output`.
  
  if (is.null(predict.vars)) {
    predict.vars <- names(output)
  }
  
  to.return <- matrix(NA, nrow=nrow(new.data), ncol=length(predict.vars))
  to.return <- as.data.frame(to.return)
  colnames(to.return) <- predict.vars
  
  for (i in 1:length(predict.vars)) {
    col <- predict.vars[i]
    scores <- NULL
    if (! col %in% colnames(new.data)) {
      next 
    }
    obj <- output[[col]]
    repeats <- length(obj$models)
    
    w.scores <- matrix(NA, nrow(new.data), ncol=repeats)
    for (j in 1:length(obj$models)) {
      model <- obj$models[[j]]
      resid <- obj$residuals[j]
      predicts.test <- predict(model, new.data)
      residuals.test <- new.data[[col]] - predicts.test 
      wscore <- residuals.test / resid
      w.scores[, j] <- wscore
    }
    
    w.scores.mean <- rowMeans(w.scores)
    to.return[, i] <- w.scores.mean
  }
  
  return (to.return)
}