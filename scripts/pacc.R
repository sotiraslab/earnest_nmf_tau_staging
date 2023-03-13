
compute.pacc <- function(df, pacc.columns,
                         cn.mask, min.required = 2,
                         higher.better = NULL) {
  cn.data <- df[cn.mask, ]
  n = nrow(df)
  k = length(pacc.columns)
  
  normed.scores <- matrix(data=NA, nrow=n, ncol=k)
  normed.scores <- as.data.frame(normed.scores)
  colnames(normed.scores) <- pacc.columns
  
  if (is.null(higher.better)) {
    higher.better <- rep(TRUE, k)
  }
  
  for (i in 1:k) {
    col <- pacc.columns[i]
    mu <- mean(cn.data[[col]], na.rm = T)
    s <- sd(cn.data[[col]], na.rm = T)
    z <- (df[[col]] - mu) / s
    
    if (! higher.better[i]) {
      z <- (-1 * z)
    }
    
    normed.scores[, i] <- z
  }
  
  pacc.score <- rowMeans(normed.scores, na.rm = T)
  count.present <- rowSums(! is.na(normed.scores))
  pacc.score <- ifelse(count.present >= min.required, pacc.score, NA)
  
  return(pacc.score)
}