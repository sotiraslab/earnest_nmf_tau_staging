sh(library(ggsignif))
sh(library(tidyverse))

my.anova <- function(x, y, data, correction='fdr', print = F) {
  fml <- as.formula(sprintf('%s ~ %s', y, x))
  anova <- aov(fml, data = data)
  posthoc <- as.data.frame(TukeyHSD(anova, method = correction)[[x]])
  posthoc.res <- posthoc %>%
    rownames_to_column('comparison') %>%
    mutate(annotation = cut(`p adj`,
                            breaks = c(0, 0.001, 0.01, 0.05, Inf),
                            labels = c('***', "**", "*", ""),
                            include.lowest = T)
    )
  if (print) {
    print(fml)
    print(summary(anova))
  }
  return (posthoc.res)
}

get_geom_sig_y_positions <- function(posthoc.res, ordered.x,
                                     start.pos = 1, gap = 1){
  
  # preprocessing
  posthoc.sig <- filter(posthoc.res, `p adj` < 0.05)
  comparisons <- str_split(posthoc.sig$comparison, '-')
  
  # data holders
  extents <- list()
  y_position <- rep(NA, nrow(posthoc.sig))
  
  # main loop
  for (i in seq_along(comparisons)){
    x <- comparisons[[i]]
    positions <- match(x, ordered.x)
    
    ypos <- NA
    check.pos <- start.pos
    while (is.na(ypos)) {
      
      # first iteration
      if (! as.character(check.pos) %in% names(extents)) {
        ypos <- check.pos
        extents[[as.character(check.pos)]] = sort(positions)
        next
      }
      
      # if entry already found, check if we can still fit new bar
      existing <- extents[[as.character(check.pos)]]
      if (max(existing) < min(positions)) {
        ypos <- check.pos
        both <- c(existing, positions)
        extents[[as.character(check.pos)]] = c(min(both), max(both))
      }
      check.pos <- check.pos + gap
    }
    
    # record found position
    y_position[i] <- ypos
  }
  
  return (y_position)
}

anova.plot <- function(x, y, data, colors, correction='fdr',
                       sig.y.start = 1, sig.y.gap = 1) {
  
  # get anova stats
  posthoc.res <- my.anova(x = x, y = y, data = data, correction = correction, print = T)
  posthoc.sig <- filter(posthoc.res, `p adj` < 0.05)
  comparisons <- str_split(posthoc.sig$comparison, '-')
  n.sig <- nrow(posthoc.sig)
  
  # get means by group
  n.categories <- length(unique(data[[x]]))
  means <- group_by(data, !!sym(x)) %>%
    summarise(Mean = mean(!!sym(y), na.rm=T)) %>%
    mutate(x=seq(0.5, by=1, length.out=n.categories),
           xend=seq(1.5, by=1, length.out=n.categories))
  
  # plot
  y_position = get_geom_sig_y_positions(posthoc.res = posthoc.res,
                                        ordered.x = sort(unique(data[[x]])),
                                        start.pos = sig.y.start,
                                        gap = sig.y.gap)
  
  ggplot(data = data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(x))) +
    geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
    geom_segment(data=means, aes(x=x, xend=xend, y=Mean, yend=Mean),
                 color='black',
                 linewidth=1) + 
    scale_fill_manual(values=colors) +
    theme_light() +
    theme(legend.position = 'none',
          text = element_text(size=20),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    coord_cartesian(ylim = c(-10, 6)) +
    scale_y_continuous(breaks=c(0, -2.5, -5, -7.5, -10)) +
    geom_signif(comparisons=comparisons,
                annotations = posthoc.sig$annotation,
                y_position = y_position,
                tip_length = 0.01,
                size=.75,
                textsize = 7) +
    xlab(x) +
    ylab(y)
}