
sh <- suppressPackageStartupMessages

sh(library(tidyverse))

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
