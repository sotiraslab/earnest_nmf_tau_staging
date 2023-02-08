
library(anticlust)

repeated.wscore <- function(control.data, test.data, y, covariates,
                            match.continuous, match.categorical,
                            portion.train = 0.8, repeats = 200, seed = T,
                            verbose = T) {

  # returned at top level
  big.output <- list()
  return.splits <- matrix(NA, nrow=nrow(control.data), ncol=repeats)
  mean.w.scores <- matrix(NA, nrow=nrow(test.data), ncol=length(y))
  
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
    w.scores <- matrix(NA, nrow=nrow(test.data), ncol=repeats)
    
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
      
      # now calculate w scores
      predicts.test <- predict(model, test.data)
      residuals.test <- test.data[[col]] - predicts.test
      wscore <- residuals.test / residual.control.sd
      
      w.scores[, i] <- wscore
      
    }
    w.scores.mean <- rowMeans(w.scores)
    collect = list(
      w.scores.mean = w.scores.mean,
      w.scores.all = w.scores,
      models = holder.models,
      residuals = holder.residuals
    )
    big.output[[col]] <- collect
    mean.w.scores[, j] <- w.scores.mean
  }
  colnames(mean.w.scores) <- y
  big.output[['..splits']] <- return.splits
  big.output[['..wscores']] <- mean.w.scores
  
  return (big.output)
}

repeated.wscore.predict <- function(output, new.data) {
  
  predict.vars <- names(output)
  predict.vars <- predict.vars[! grepl('^\\.\\.', predict.vars, perl = T)]
  
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