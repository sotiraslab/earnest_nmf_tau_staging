# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(dplyr))
sh(library(ggplot2))
sh(library(reshape))
sh(library(tidyr))

stage.heatmap <- function(df, cols=NULL, return.data=F, omit.zero=F) {
  
  if (is.null(cols)) {
    cols <- colnames(df)
  }
  
  df <- df[, cols]
  ordered.by.sum <- names(sort(colSums(df), decreasing = T))
  binarized.cols.df <- df[, c(ordered.by.sum)]
  binarized.cols.df <- binarized.cols.df[do.call(order, binarized.cols.df ), ]
  if (omit.zero) {
    binarized.cols.df <- binarized.cols.df[rowSums(binarized.cols.df[cols]) > 0, ]
  }
  idx <- 1:nrow(binarized.cols.df)
  binarized.cols.df$Idx <- idx
  
  if (return.data) return (binarized.cols.df)
  
  totals <- colSums(binarized.cols.df[, cols])
  totals <- data.frame(variable = names(totals), value = unname(totals),
                       Idx=0)
  
  hmap.data <- as.data.frame(binarized.cols.df) %>%
    reshape::melt(id.vars = 'Idx', measure.vars = cols) %>%
    mutate(variable = factor(variable, levels=ordered.by.sum))
  
  p <- ggplot(hmap.data, aes(x=variable, y=Idx, fill=value)) +
    geom_tile(color='black') +
    scale_fill_gradient(low = "white", high = "red", limits=c(0,1)) +
    geom_text(data=totals, size=3, nudge_y=5, aes(label=value)) +
    theme_classic() +
    theme(legend.position = 'none', line=element_blank(),
          axis.text.x = element_text(size=6, angle=45, vjust=1.15, hjust=1),
          axis.ticks.length = unit(0, 'cm')) +
    scale_y_reverse() +
    ylab('Participant') +
    xlab('Region')
  
  return (p)
}

stage.heatmap.by <- function(df, by, cols=NULL, cats=NULL, return.data=F, omit.zero=F,
                             colors=NULL, empty.in.legend=T, empty.name='Negative') {
  
  # set defaults
  if (is.null(cols)) {
    cols <- colnames(df)
  }
  if (is.null(cats)) {
    cats <- as.character(sort(unique(df[[by]])))
  }
  
  # trim df
  all.cols <- c(by, cols)
  df <- df[, all.cols]
  
  # get order of columns, sort columns by that order
  # place by column in front
  ordered.by.sum <- names(sort(colSums(df[, cols]), decreasing = T))
  binarized.cols.df <- df[, c(by, ordered.by.sum)]
  
  # limit data to only include cats
  binarized.cols.df <- binarized.cols.df[binarized.cols.df[[by]] %in% cats, ]
  
  # collect totals
  totals <- colSums(binarized.cols.df[, cols])
  totals <- data.frame(variable = names(totals), value = unname(totals),
                       Idx=0)
  
  # scale cells with the by column
  binarized.cols.df[[by]] <- factor(binarized.cols.df[[by]], levels=cats)
  binarized.cols.df[, cols] <- as.numeric(binarized.cols.df[[by]]) * binarized.cols.df[, cols]
  groupsizes <- count(binarized.cols.df, across(by))
  ncats <- nrow(groupsizes)
  
  # sort rows with by column, then positivity
  binarized.cols.df <- binarized.cols.df[do.call(order, binarized.cols.df), ]
  
  # omit zeros if request
  if (omit.zero) {
    binarized.cols.df <- binarized.cols.df[rowSums(binarized.cols.df[cols]) > 0, ]
  }
  idx <- 1:nrow(binarized.cols.df)
  binarized.cols.df$Idx <- idx
  
  if (return.data) return (binarized.cols.df)
  
  t <- hmap.data <- as.data.frame(binarized.cols.df) %>%
    reshape::melt(id.vars = 'Idx', measure.vars = cols)
  
  hmap.data <- as.data.frame(binarized.cols.df) %>%
    reshape::melt(id.vars = 'Idx', measure.vars = cols) %>%
    mutate(variable = factor(variable, levels=ordered.by.sum), !!by := factor(value, labels=c(empty.name, cats)))
  
  # if (is.null(colors)) {
  #   colors = c('0'='white')
  #   v <- rocket(ncats)
  #   for (i in 1:ncats) {
  #     colors[[as.character(i)]] <- v[[i]]
  #   }
  # } else {
  #   colors = 
  # }
  
  plotcolors <- list()
  plotcolors[[ empty.name ]] <- 'white'
  if (is.null(colors)) {
    colors <- colormap('viridis', nshades = ncats)
  } 
  for (i in 1:ncats) {
    plotcolors[[cats[i]]] <- colors[[i]]
  }
  
  if (empty.in.legend) breaks <- c(empty.name, cats) else breaks <- cats 
  
  p <- ggplot() +
    geom_tile(data =hmap.data,  color='black', aes_string(x='variable', y='Idx', fill=by)) +
    scale_fill_manual(values=plotcolors, breaks=breaks) +
    geom_text(data=totals, size=3, nudge_y=nrow(df)/80, aes(x=variable, y=Idx, label=value)) +
    theme_classic() +
    theme(line=element_blank(),
          axis.text.x = element_text(size=6, angle=45, vjust=1.15, hjust=1),
          axis.ticks.length = unit(0, 'cm')) +
    scale_y_reverse() +
    ylab('Participant') +
    xlab('Region')
  
  ticker <- 0
  for (i in 1:(ncats-1)) {
    x <- groupsizes$n[i]
    ticker <- ticker + as.numeric(groupsizes$n[i])
    p <- p + geom_hline(yintercept = ticker + 0.5, color='red')
  }
  
  
  return (p)
}

stage.results <- function(df) {
  df <- df[, names(sort(colSums(df), decreasing = T))]
  n <- nrow(df)
  
  stagedf <- data.frame(Region=colnames(df),
                       NPos=colSums(df))
  rownames(stagedf) <- NULL
  
  row.sums <- rowSums(df)
  nopathology <- df[row.sums == 0, ]
  pathology <- df[row.sums > 0, ]
  
  diffs <- apply(pathology, 1, diff)
  instage <- apply(diffs <= 0, 2, all)
  
  results = list(stagedf = stagedf, 
                 not.positive = nrow(nopathology),
                 any.positive = nrow(pathology),
                 in.stage = sum(instage), 
                 out.stage = sum(! instage))
  
  return (results)
}
