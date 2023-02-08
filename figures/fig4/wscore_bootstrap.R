# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(gtools))
sh(library(ggplot2))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =======

PATH.ADNI.WSCORES <- '../fig3/adni_data_with_wscores.csv'
PATH.ADNI.ORDER <- '../fig3/adni_wscore_stage_order.csv'

# === Read data ========

df <- read.csv(PATH.ADNI.WSCORES)

df <- df %>%
  filter(Group == 'TrainingBaseline')

# === Get W scored region order =========

w.stages <- read.csv(PATH.ADNI.ORDER)
w.order <- w.stages$Region
w.cols <- paste(w.order, '.WScore', sep='')

nice.names <- gsub('Cmp.|.WScore', '', w.cols)

df.w <- df[, w.cols]
colnames(df.w) <- nice.names
df.w <- as.data.frame(ifelse(df.w >= 2.5, 1, 0))

# === List all comparisons ========

comparisons <- t(combn(nice.names, 2))
colnames(comparisons) <- c('A', 'B')
n.compare <- nrow(comparisons)

# === Calculate observed statistics

observed.pos <- colSums(df.w)

observed.diffs <- as.data.frame(comparisons)
observed.diffs$Diff <- NA
for (i in 1:nrow(observed.diffs)) {
  a <- observed.diffs[i, 'A']
  b <- observed.diffs[i, 'B']
  observed.diffs[i, 'Diff'] <- observed.pos[a] - observed.pos[b]
}

# === Bootstrapping =======

set.seed(42)
N <- 5000

# Create holder for nulls for all comparisons
nulls <- matrix(data = NA, nrow = n.compare, ncol = N)

# run
for (n in 1:N) {
  idx <- sample(1:nrow(df.w), size = nrow(df.w), replace = T)
  boot <- df.w[idx, ]
  boot.pos <- colSums(boot)
  null.pos <- boot.pos - observed.pos
  
  for (i in 1:nrow(comparisons)) {
    a <- comparisons[i, 'A']
    b <- comparisons[i, 'B']
    nulls[i, n] <- null.pos[a] - null.pos[b]
  }
}

# === P-values ========

observed.diffs$p <- rowMeans(nulls >= observed.diffs$Diff)
observed.diffs$p.corrected <- p.adjust(observed.diffs$p, method = 'fdr')
observed.diffs$log.p.corrected <- -log10(observed.diffs$p.corrected)

# === Collect all pairwise comparisons ========

all.pairs <- as.data.frame(permutations(length(nice.names), 2, nice.names, repeats.allowed = T))
colnames(all.pairs) <- c('A', 'B')

all.pairs <- left_join(all.pairs, observed.diffs, by=c('A', 'B'))
all.pairs$Main <- ! is.na(all.pairs$log.p.corrected)

# ==== Plot  =========

plot.data <- all.pairs %>%
  mutate(plot.p = ifelse(log.p.corrected >= 5, 5, log.p.corrected),
         plot.p = ifelse(plot.p<= -log10(0.05), 0, plot.p),
         plot.p = ifelse(is.na(plot.p), 0, plot.p),
         A = factor(A, levels=rev(colnames(df.w))),
         B = factor(B, levels=colnames(df.w)))

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
  ylab('PTC') +
  xlab('PTC')

ggsave('adni_wscore_bootstrap.png', width = 8, height = 8)

# === Plot showing stages =========

plot.data <- all.pairs %>%
  mutate(plot.p = ifelse(log.p.corrected >= 5, 5, log.p.corrected),
         plot.p = ifelse(plot.p<= -log10(0.05), 0, plot.p),
         plot.p = ifelse(is.na(plot.p), 0, plot.p),
         A = factor(A, levels=rev(colnames(df.w))),
         B = factor(B, levels=colnames(df.w)))

# this was for creating rectangles over the lower triangle to indicate stages
# seemed a little too messy

# stage.rects <- data.frame(xmin = c(1, 2, 5, 7) - 0.5,
#                           xmax = c(1, 4, 6, 8) + 0.5,
#                           ymin = c(1, 1, 1, 1) - 0.5,
#                           ymax = c(8, 7, 4, 2) + 0.5,
#                           Stage = factor(1:4))

stage.geoms <- data.frame(xmin = c(1, 2, 5, 7) - 0.4,
                          xmax = c(1, 4, 6, 8) + 0.4,
                          ymin = rep(8.75, 4),
                          ymax = rep(8.9, 4),
                          Stage = factor(1:4),
                          names = paste('Stage', 1:4),
                          xcenter = c(1, 3, 5.5, 7.5),
                          ycenter = rep(9.1, 4))

# https://mikemol.github.io/technique/colorblind/2018/02/11/color-safe-palette.html
stage.colors <- c(`1` = '#009E73', `2` = '#F0E442', `3` = '#E69F00', `4` = '#D55E00')

ggplot() +
  geom_tile(data = plot.data, aes(x = B, y = A, fill = plot.p), linewidth=1.5, color='white') +
  scale_fill_colormap(colormap='viridis') +
  coord_equal() +
  labs(fill="-log10(p)") +
  ggnewscale::new_scale_fill() + 
  geom_rect(data = stage.geoms,
            mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Stage),
            color='black') +
  geom_text(data = stage.geoms, aes(x = xcenter, y = ycenter, label = names),
            angle=45,
            hjust=0,
            fontface='bold') + 
  scale_fill_manual(values = stage.colors, guide='none') +
  scale_x_discrete(expand=expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand=expansion(mult = c(0, 0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle=30, hjust=1),
        text = element_text(size=15)) +
  ylab('PTC') +
  xlab('PTC')

ggsave('adni_wscore_bootstrap_stage_labeled.png', width = 8, height = 8)