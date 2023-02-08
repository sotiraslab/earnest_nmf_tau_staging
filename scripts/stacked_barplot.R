
sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(tidyverse))

stacked.barplot <- function(df, xcol, ycol, levels=NULL, colors=NULL) {
  
  # create data
  plot.data <- df %>%
    group_by_at(c(xcol, ycol)) %>%
    summarise(N=n())
  
  if (! is.null(levels)) {
    plot.data <- filter(plot.data, !!sym(xcol) %in% levels)
  }
  plot.data <- plot.data %>%
    ungroup() %>%
    group_by(!!sym(xcol)) %>%
    mutate(Percent = N / sum(N) * 100)
  
  # get totals for text
  group.sums <- plot.data %>%
    group_by(!!sym(xcol)) %>%
    summarise(total=sum(N)) %>%
    ungroup()
  
  # plot
  p <- ggplot() +
    geom_bar(data = plot.data, aes(fill=!!sym(ycol), y=Percent, x=!!sym(xcol)), stat="identity", color='black') +
    geom_text(data = group.sums, aes(x=!!sym(xcol), y=105, label=total), size=6) +
    theme_classic() +
    ylab('Observations (%)') +
    xlab(xcol) +
    scale_y_continuous(expand=expansion(mult=c(0, .1)), breaks = c(0, 25, 50, 75, 100)) +
    guides(fill=guide_legend(title="Stage")) +
    theme(text = element_text(size=20),
          axis.line.y = element_blank()) +
    geom_segment(aes(y=0,yend=100,x=-Inf,xend=-Inf), color='black', linewidth=1)
  
  if (! is.null(colors)) {
    p <- p + scale_fill_manual(values=colors)
  }
  
  return(p)
}