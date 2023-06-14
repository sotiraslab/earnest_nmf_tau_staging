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

PATH.OASIS.WSCORES <- '../../derivatives/oasis3/data_with_wscores.csv'
PATH.OASIS.ORDER <- '../../derivatives/oasis3/wscore_stage_order.csv'

# === Read data ========

df <- read.csv(PATH.OASIS.WSCORES)

df <- df %>%
  filter(Group == 'TrainingSet')

# === Get W scored region order =========

w.stages <- read.csv(PATH.OASIS.ORDER)
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

ggsave('oasis3_wscore_bootstrap.png', width = 8, height = 8)
