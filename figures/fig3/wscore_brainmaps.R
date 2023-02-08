# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggplot2))
sh(library(ggseg))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# ==== Required files ========

PATH.ADNI.WSCORES <- 'adni_data_with_wscores.csv'
PATH.ADNI.ORDER <- 'adni_wscore_stage_order.csv'
PATH.OASIS.WSCORES <- 'oasis_data_with_wscores.csv'
PATH.OASIS.ORDER <- 'oasis_wscore_stage_order.csv'
PATH.REGION.ASSIGNMENT <- 'assignments.csv'

# ==== Plot adni ====

adni.pos <- read.csv(PATH.ADNI.ORDER)
assignments <- read.csv(PATH.REGION.ASSIGNMENT)

adni.w <- read.csv(PATH.ADNI.WSCORES) %>%
  filter(Group == 'TrainingBaseline') %>%
  select(contains('WScore')) %>%
  mutate(across(.fns = function (x) ifelse(x >= 2.5, 1, 0)))

any.pos <- sum(rowSums(adni.w) > 0)

merger <- adni.pos %>%
  select(Region, NPos) %>%
  mutate(NPos = NPos / any.pos * 100)
colnames(merger) <- c('name', 'npos')

plot.data <- left_join(assignments, merger, by='name')

p <- ggplot(plot.data) +
  geom_brain(atlas = dk,
             aes(fill=npos),
             position = position_brain(hemi ~ side )) +
  theme_void() + 
  scale_fill_gradient(low='white', high='red', limits=c(20, 90)) +
  labs(fill = '% Elevated') +
  ggtitle('ADNI-ADS')
p

ggsave('adni_wscore_map.png', p, width=7, height=7)

write.csv(plot.data, 'adni_wscore_map.csv', row.names = F)

# ---- Plot ADRC -------

oasis.pos <- read.csv(PATH.OASIS.ORDER)

oasis.w <- read.csv(PATH.OASIS.WSCORES) %>%
  filter(Group == 'TrainingSet') %>%
  select(contains('WScore')) %>%
  mutate(across(.fns = function (x) ifelse(x >= 2.5, 1, 0)))

any.pos <- sum(rowSums(oasis.w) > 0)

merger <- oasis.pos %>%
  select(Region, NPos) %>%
  mutate(NPos = NPos / any.pos * 100)
colnames(merger) <- c('name', 'npos')
merger$name <- gsub('.WScore', '', merger$name)

plot.data <- left_join(assignments, merger, by='name')

p <- ggplot(plot.data) +
  geom_brain(atlas = dk,
             aes(fill=npos),
             position = position_brain(hemi ~ side )) +
  theme_void() + 
  scale_fill_gradient(low='white', high='red', limits=c(20, 90)) +
  labs(fill = '% Elevated') +
  ggtitle('OASIS3-ADS')
p

ggsave('oasis_wscore_map.png', p, width=7, height=7)

write.csv(plot.data, 'oasis_wscore_map_data.csv', row.names = F)

