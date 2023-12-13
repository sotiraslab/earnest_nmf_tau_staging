# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggsignif))
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

# ==== plot 1: positivity by cutoff ========

regional_positivity <- function(data, cutoff) {
  thr <- ifelse(data >= cutoff, 1, 0)
  return (colSums(thr))
}

values <- sapply(cutoffs, function(x) regional_positivity(wscores, x))
colnames(values) <- cutoffs
plot.data <- values %>%
  as.data.frame() %>%
  rownames_to_column('PTC') %>%
  pivot_longer(-PTC, names_to = 'threshold', values_to = 'npos') %>%
  mutate(PTC = factor(PTC, levels=ptc.order),
         threshold = as.numeric(threshold))

ggplot(data=plot.data, aes(x=threshold, y=npos, group=PTC, color=PTC, linetype=PTC)) +
  geom_line(linewidth=1) +
  geom_point() + 
  theme_light() +
  ylab('Tau-positive') +
  xlab('Cutoff (W-score)') +
  theme(text = element_text(size=15),
        legend.position = 'bottom') +
  guides(color = guide_legend(nrow=3))
  

ggsave('wscore_order_lineplot.png', width = 8, height = 6, units = 'in')

# ==== plot 2: heatmap of ordering ========

regional_ranks <- function(data, cutoff) {
  thr <- ifelse(data >= cutoff, 1, 0)
  npos <- colSums(thr)
  ranks <- rank(-npos, ties.method = 'min')
  return (ranks)
}

values <- sapply(cutoffs, function(x) regional_ranks(wscores, x))
colnames(values) <- cutoffs
plot.data <- values %>%
  as.data.frame() %>%
  rownames_to_column('PTC') %>%
  pivot_longer(-PTC, names_to = 'threshold', values_to = 'Rank') %>%
  mutate(PTC = factor(PTC, levels=rev(ptc.order)),
         Rank = factor(Rank))

colors = colormap('YIOrRd', nshades = 8)

ggplot(data = plot.data, aes(x = threshold, y = PTC,
                             fill = Rank)) +
  geom_tile(color='black') +
  theme_light() +
  coord_equal() +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0),
                   breaks = seq(2, 3, .1)) +
  scale_fill_manual(values=colors) +
  xlab('Cutoff (W-score)')

ggsave('wscore_order_heatmap.png', width = 8, height = 3, units = 'in')