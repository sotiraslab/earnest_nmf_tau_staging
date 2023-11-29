# === Import =======

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggseg))
sh(library(R.matlab))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set working directory =======

setwd(this.dir())

# === Required Files =======

PATH.WTA <- '../../derivatives/adni/ptc8_winner_take_all.csv'
PATH.PTC.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'
PATH.NEWSTAGES <- '../../derivatives/adni/ptc8_stage_assignment.csv'

# ==== read =======

wta <- read.csv(PATH.WTA)
ptc.order <- read.csv(PATH.PTC.ORDER)
ptc.stages <- read.csv(PATH.NEWSTAGES)

# === make binary PTC assignment df ========

ptc.mat <- matrix(NA, nrow = nrow(wta), ncol = nrow(ptc.order))

for (i in seq_along(ptc.order$Region)) {
  ptc <- ptc.order$Region[i]
  ptc.mat[, i] <- wta$name == ptc
}

rownames(ptc.mat) <- wta$label
colnames(ptc.mat) <- ptc.order$Region

# === make binary PTCStaging assignment df =======

stage.mat <- matrix(NA, nrow = nrow(wta), ncol = 4)

for (i in 1:4) {
  stage.mat[, i] <- ptc.stages$value == i
}

rownames(stage.mat) <- wta$label
colnames(stage.mat) <- 1:4

# ==== make braak assignments =======

braak1.regs <- c('ENTORHINAL')

braak3.regs <- c('PARAHIPPOCAMPAL',
                 'FUSIFORM',
                 'LINGUAL',
                 'AMYGDALA')

braak4.regs <- c('MIDDLETEMPORAL',
                 'CAUDALANTERIORCINGULATE',
                 'ROSTRALANTERIORCINGULATE',
                 'POSTERIORCINGULATE',
                 'ISTHMUSCINGULATE',
                 'INSULA',
                 'INFERIORTEMPORAL',
                 'TEMPORALPOLE')

braak5.regs <- c('SUPERIORFRONTAL',
                 'LATERALORBITOFRONTAL',
                 'MEDIALORBITOFRONTAL',
                 'FRONTALPOLE',
                 'CAUDALMIDDLEFRONTAL',
                 'ROSTRALMIDDLEFRONTAL',
                 'PARSOPERCULARIS',
                 'PARSORBITALIS',
                 'PARSTRIANGULARIS',
                 'LATERALOCCIPITAL',
                 'SUPRAMARGINAL',
                 'INFERIORPARIETAL',
                 'SUPERIORTEMPORAL',
                 'SUPERIORPARIETAL',
                 'PRECUNEUS',
                 'BANKSSTS',
                 'TRANSVERSETEMPORAL')

braak6.regs <- c('PERICALCARINE',
                 'POSTCENTRAL',
                 'CUNEUS',
                 'PRECENTRAL',
                 'PARACENTRAL')

braak.mat <- matrix(NA, nrow = nrow(wta), ncol = 5)

braaks <- list(braak1.regs,
               braak3.regs,
               braak4.regs,
               braak5.regs,
               braak6.regs)

for (i in 1:5) {
  regions <- braaks[[i]]
  pattern <- paste(tolower(regions), collapse='|')
  braak.mat[, i] <- str_detect(wta$label, pattern)
}

rownames(braak.mat) <- wta$label
colnames(braak.mat) <- c('BraakI', 'BraakIII', 'BraakIV', 'BraakV', 'BraakVI')

# === Similarity functions =======

# https://www.r-bloggers.com/2021/11/how-to-calculate-jaccard-similarity-in-r/
jaccard <- function(a, b) {
  intersection = sum(a & b)
  union = sum(a | b)
  return (intersection/union)
}

dice <- function(a, b) {
  top <- 2 * sum(a & b)
  bottom <- sum(a) + sum(b)
  return (top/bottom)
}

# === Similarity for PTCs and Braak ROIs =======

sim.mat <- matrix(data = NA, nrow = ncol(braak.mat), ncol = ncol(ptc.mat))

for (i in 1:ncol(braak.mat)) {
  for (j in 1:ncol(ptc.mat))
    sim.mat[i, j] <- dice(braak.mat[, i], ptc.mat[, j])
}

colnames(sim.mat) <- colnames(ptc.mat)
rownames(sim.mat) <- colnames(braak.mat)

plot.data <- sim.mat %>%
  as.data.frame() %>%
  rownames_to_column('Braak') %>%
  pivot_longer(-Braak, names_to = 'PTC', values_to = 'Dice') %>%
  mutate(Braak = str_replace(Braak, 'Braak',''),
         PTC = str_replace(PTC, 'Cmp.', ''),
         Braak = factor(Braak, levels=c('VI', 'V', 'IV', 'III', 'I')),
         PTC = factor(PTC, levels = str_replace(ptc.order$Region, 'Cmp.', '')))

plot <- ggplot(plot.data, aes(x = PTC, y = Braak, fill = Dice)) + 
  geom_tile() + 
  coord_equal() +
  scale_x_discrete(expand=expansion(mult = c(0, 0))) +
  scale_y_discrete(expand=expansion(mult = c(0, 0))) + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size=12)) + 
  scale_fill_colormap(colormap='velocity-blue', limits=c(0, 1))

ggsave(plot = plot, filename = 'braak_ptc_similarity.png', width=8, height=4)

plot

# === Similarity for PTC-staging and Braak ROIs =======

sim.mat <- matrix(data = NA, nrow = ncol(braak.mat), ncol = ncol(stage.mat))

for (i in 1:ncol(braak.mat)) {
  for (j in 1:ncol(stage.mat))
    sim.mat[i, j] <- dice(braak.mat[, i], stage.mat[, j])
}

colnames(sim.mat) <- colnames(stage.mat)
rownames(sim.mat) <- colnames(braak.mat)

plot.data <- sim.mat %>%
  as.data.frame() %>%
  rownames_to_column('Braak') %>%
  pivot_longer(-Braak, names_to = 'PTCStage', values_to = 'Dice') %>%
  mutate(Braak = str_replace(Braak, 'Braak',''),
         Braak = factor(Braak, levels=c('VI', 'V', 'IV', 'III', 'I')))

plot <- ggplot(plot.data, aes(x = PTCStage, y = Braak, fill = Dice)) + 
  geom_tile() + 
  coord_equal() +
  scale_x_discrete(expand=expansion(mult = c(0, 0))) +
  scale_y_discrete(expand=expansion(mult = c(0, 0))) + 
  theme(text = element_text(size=12)) + 
  scale_fill_colormap(colormap='velocity-blue', limits=c(0, 1))

ggsave(plot = plot, filename = 'braak_newstaging_similarity.png', width=8, height=4)

plot