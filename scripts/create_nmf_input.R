create.nmf.input <- function(df, features, outdata, outids, outfeatures=NULL,
                             id.cols=c('Subject')) {
  identifiers <- df[, id.cols]
  data <- as.matrix(df[features])
  X <- t(data)
  
  cat("NAs:", sum(is.na(X)), '\n')
  cat("Negatives:", sum(X < 0), '\n')
  
  write.table(X, outdata, row.names = F, col.names = F, sep=',')
  write.csv(identifiers, outids, quote = F)
  if (! is.null(outfeatures)) {
    feature.df <- data.frame(Feature=features)
    write.csv(feature.df, outfeatures, quote = F)
  }
  
  return(X)
}