# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggsignif))
sh(library(gtools))
sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.DATA <- '../../derivatives/adni/data_with_wscores.csv'
PATH.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'

# === read  ========

df.all <- read.csv(PATH.DATA)

# select only training data
df <- df.all %>%
  filter(Group == 'TrainingBaseline')

wscores <- df %>%
  select(contains('WScore'))
colnames(wscores) <- str_replace_all(colnames(wscores), 'Cmp.|.WScore', '')

ptc.order <- read.csv(PATH.ORDER)$Region
ptc.order <- str_replace_all(ptc.order, 'Cmp.|.WScore', '')

# === cutoffs to assess =========

cutoffs <- seq(2.0, 3.0, .05)

# === create function for bootstrap staging =======

bootstrap.staging <- function(binary.df, N=5000, seed=NULL) {
  observed.pos <- colSums(binary.df)
  order.regions <- names(sort(observed.pos, decreasing = T))
  
  comparisons <- t(combn(order.regions, 2))
  colnames(comparisons) <- c('A', 'B')
  n.compare <- nrow(comparisons)
  
  
  
  observed.diffs <- as.data.frame(comparisons)
  observed.diffs$Diff <- NA
  for (i in 1:nrow(observed.diffs)) {
    a <- observed.diffs[i, 'A']
    b <- observed.diffs[i, 'B']
    observed.diffs[i, 'Diff'] <- observed.pos[a] - observed.pos[b]
  }
  
  if (! is.null(seed)) {
    set.seed(42)
  }
  
  # Create holder for nulls for all comparisons
  nulls <- matrix(data = NA, nrow = n.compare, ncol = N)
  
  # run
  for (n in 1:N) {
    idx <- sample(1:nrow(binary.df), size = nrow(binary.df), replace = T)
    boot <- binary.df[idx, ]
    boot.pos <- colSums(boot)
    null.pos <- boot.pos - observed.pos
    
    for (i in 1:nrow(comparisons)) {
      a <- comparisons[i, 'A']
      b <- comparisons[i, 'B']
      nulls[i, n] <- null.pos[a] - null.pos[b]
    }
  }
  
  observed.diffs$p <- rowMeans(nulls >= observed.diffs$Diff)
  observed.diffs$p.corrected <- p.adjust(observed.diffs$p, method = 'fdr')
  observed.diffs$log.p.corrected <- -log10(observed.diffs$p.corrected)
  
  return (observed.diffs)
}

infer.stages <- function(output) {
  regions <- unique(c(output$A, output$B))
  n.regions <- length(regions)
  
  all.pairs <- as.data.frame(permutations(length(regions), 2, regions, repeats.allowed = T))
  colnames(all.pairs) <- c('A', 'B')
  
  all.pairs <- left_join(all.pairs, output, by=c('A', 'B'))
  all.pairs$Main <- ! is.na(all.pairs$log.p.corrected)
  
  mat.data <- all.pairs %>%
    mutate(A = factor(A, levels=regions),
           B = factor(B, levels=regions),
           significant = ifelse(p.corrected < 0.05, 1, 0),
           significant = ifelse(is.na(significant), 0, significant)) %>%
    select(A, B, significant) %>% 
    arrange(A, B) %>%
    pivot_wider(id_cols = A, names_from = B, values_from = significant) %>%
    select(-A) %>%
    as.matrix()
  
  rowdiff <- apply(mat.data, 1, diff)
  coldiff <- apply(mat.data, 2, diff)
  if (any(rowdiff < 0) | any(coldiff > 0)) {
    warning('Cannot infer stages automatically.  Retuning NA')
    return (rep(NA, n.regions))
  }
  
  
  stages <- as.integer(factor(colSums(mat.data)))
  return (stages)
}

bootstrap.stage.heatmap <- function(output) {
  
  regions <- unique(c(output$A, output$B))
  
  all.pairs <- as.data.frame(permutations(length(regions), 2, regions, repeats.allowed = T))
  colnames(all.pairs) <- c('A', 'B')
  
  all.pairs <- left_join(all.pairs, output, by=c('A', 'B'))
  all.pairs$Main <- ! is.na(all.pairs$log.p.corrected)
  
  # ==== Plot  =========
  
  plot.data <- all.pairs %>%
    mutate(plot.p = ifelse(log.p.corrected >= 5, 5, log.p.corrected),
           plot.p = ifelse(plot.p <= -log10(0.05), 0, plot.p),
           plot.p = ifelse(is.na(plot.p), 0, plot.p),
           A = factor(A, levels=rev(regions)),
           B = factor(B, levels=regions))
  
  ggplot(data = plot.data, aes(x = B, y = A, fill = plot.p)) +
    geom_tile(linewidth=1.5, color='white') +
    coord_equal() +
    scale_fill_colormap(colormap='viridis') +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#440154ff"),
          axis.text.x = element_text(angle=30, hjust=1),
          text = element_text(size=15)) +
    labs(fill="-log10(p)") + 
    ylab('RegionA') +
    xlab('RegionB')
}

pipeline <- function(wscores, cutoff) {
  print(cutoff)
  binary.df <- as.data.frame(ifelse(wscores >= cutoff, 1, 0))
  output <- bootstrap.staging(binary.df, N=5000, seed=42)
  stages <- infer.stages(output)
}

# ------ run staging for all thresholds -------

values <- sapply(cutoffs, function(x) pipeline(wscores, x))

# ------ plot ---------

# heatmap data
colnames(values) <- cutoffs
plot.data <- values %>%
  as.data.frame() %>%
  mutate(PTC = ptc.order) %>%
  pivot_longer(-PTC, names_to = 'threshold', values_to = 'Stage') %>%
  mutate(PTC = factor(PTC, levels=rev(ptc.order)),
         Stage = factor(Stage))
# draw
colors = c('#009E73', '#F0E442', '#E69F00', '#D55E00', '#CC79A7')

ggplot(data = plot.data, aes(x = threshold, y = PTC,
                             fill = Stage)) +
  geom_tile(color='black') +
  theme_light() +
  coord_equal() +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0),
                   breaks = seq(2, 3, .1)) +
  scale_fill_manual(values=colors) +
  xlab('Cutoff (W-score)')

ggsave('staging_sensitivity.png', width = 8, height = 3, units = 'in')

# ------- alternative plot ---------------

# draws boxes around stages

# # heatmap data
# colnames(values) <- cutoffs
# plot.data <- values %>%
#   as.data.frame() %>%
#   mutate(PTC = ptc.order) %>%
#   pivot_longer(-PTC, names_to = 'threshold', values_to = 'Stage') %>%
#   mutate(PTC = factor(PTC, levels=rev(ptc.order)),
#          Stage = factor(Stage))
# 
# # data for boxes around stages
# total.boxes <- sum(apply(values, 2, max))
# box.data <- data.frame(xmin = rep(NA, total.boxes),
#                        xmax = rep(NA, total.boxes),
#                        ymin = rep(NA, total.boxes),
#                        ymax = rep(NA, total.boxes))
# 
# row <- 1
# for (i in 1:ncol(values)) {
#   stages <- values[, i]
#   n.stages <- max(stages)
#   box.data[row:(row+n.stages-1), 1] = i - 0.5
#   box.data[row:(row+n.stages-1), 2] = i + 0.5
#   box.data[row:(row+n.stages-1), 3] = 0.5 + c(0, cumsum(rev(table(stages)))[1:(n.stages-1)])
#   box.data[row:(row+n.stages-1), 4] = 0.5 + cumsum(rev(table(stages)))
#   row <- row + n.stages
# }
# 
# # draw
# colors = c('#009E73', '#F0E442', '#E69F00', '#D55E00', '#CC79A7')
# 
# ggplot(data = plot.data, aes(x = threshold, y = PTC,
#                              fill = Stage)) +
#   geom_tile(color=NA) +
#   theme_light() +
#   coord_equal() +
#   scale_y_discrete(expand = expansion(add=.75)) +
#   scale_x_discrete(expand = expansion(add=.75),
#                    breaks = seq(2, 3, .1)) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank()) +
#   scale_fill_manual(values=colors) +
#   xlab('Cutoff (W-score)') +
#   geom_rect(data = box.data,
#             aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
#             inherit.aes = F,
#             fill=NA,
#             color='black',
#             linewidth = 1.5)
