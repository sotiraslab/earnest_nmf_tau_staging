# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.DATA <- '../../derivatives/oasis3//data_with_staging.csv'
PATH.PTC.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'
PATH.WTA <- '../../derivatives/adni/ptc8_winner_take_all.csv'
PATH.PLOTSCRIPT <- '../../scripts/ggseg_plots.R'

# === Read ========

df <- read.csv(PATH.DATA)
ptc.order <- read.csv(PATH.PTC.ORDER)$Region

df <- df %>%
  filter(Group == 'TrainingSet') %>%
  mutate(Stageable = ifelse(PTCStage %in% c('1', '2', '3', '4'), 'Stageable', NA),
         Stageable = ifelse(PTCStage == '0', 'Stage 0', Stageable),
         Stageable = ifelse(PTCStage == 'NS', 'NS', Stageable),
         Laterality = abs(Cmp.LeftParietalTemporal - Cmp.RightParietalTemporal),
         Stageable = factor(Stageable, levels = c('Stage 0', 'Stageable', 'NS')))

# === Create functions ========

lm.statistics <- function(dependent, adjust.tau = F) {
  mask.stageable <- df$Stageable == 'Stageable'
  mask.ns <- df$Stageable == 'NS'
  m.data <- df[mask.ns | mask.stageable, ]

  if (adjust.tau) {
    fml <- as.formula(sprintf('%s ~ Stageable + TotalCtxTauMean', dependent))
  } else {
    fml <- as.formula(sprintf('%s ~ Stageable', dependent))
  }

  m <- lm(fml, data = m.data)

  data.stageable <- df[mask.stageable, dependent]
  data.ns <- df[mask.ns, dependent]

  vals <- list(dependent,
            mean(data.stageable, na.rm = T),
            sd(data.stageable, na.rm = T),
            mean(data.ns, na.rm = T),
            sd(data.ns, na.rm = T),
            summary(m)$coefficients[2, 4])
}

create.stats.table <- function() {
  dependents <- c('Age',
                  'PACC.Original',
                  'MMSE',
                  'Centiloid',
                  'TotalCtxTauMean',
                  ptc.order,
                  'Laterality')
  
  mat <- matrix(NA, nrow = length(dependents), ncol = 6)
  stat.data <- as.data.frame(mat)
  colnames(stat.data) <- c('Variable',
                        'Stageable.mean',
                        'Stageable.sd',
                        'NS.mean',
                        'NS.sd',
                        'pval')
  
  for (i in seq_along(dependents)) {
    dependent <- dependents[i]
    if (dependent %in% c(ptc.order, 'Laterality')) adjust <- T else adjust <- F
    stat.data[i, ] <- lm.statistics(dependent, adjust.tau = adjust)
  }
  stat.data$pval.adjust <- p.adjust(stat.data$pval, method = 'fdr')
  return (stat.data)
}

stat.data <- create.stats.table()

# ==== Save tables ========

# raw
write.csv(stat.data, 'oasis_stats_raw.csv', quote = F, na = '', row.names = F)

# string version
str.table <- stat.data %>%
  mutate(Stageable = sprintf('%s (%s)', round(Stageable.mean, 2), round(Stageable.sd, 2)),
         NonStageable = sprintf('%s (%s)', round(NS.mean, 2), round(NS.sd, 2)),
         annotation = cut(
           pval.adjust,
           breaks = c(0, 0.001, 0.01, 0.05, Inf),
           labels = c('***', "**", "*", ""),
           include.lowest = T),
         pvalue = str_c(as.character(round(pval.adjust, 3)), annotation, sep = ' '),
         Variable = str_replace(Variable, 'Cmp.', '')) %>%
  select(Variable, Stageable, NonStageable, pvalue)

write.csv(str.table, 'oasis_stats_formatted.csv', quote = F, na = '', row.names = F)

# === Plot proportion positive ==========

source(PATH.PLOTSCRIPT)
wta <- read.csv(PATH.WTA)

w.cols <- colnames(df)[str_detect(colnames(df), '.WScore')]
w.df <- df[df$PTCStage == 'NS', w.cols]
pos.df <- as.data.frame(w.df >= 2.5)
ptc.pos <- colSums(pos.df)
merger <- data.frame(name = names(ptc.pos), npos = unname(ptc.pos))
merger$name <- str_replace(merger$name, '.WScore', '')
plot.data <- left_join(wta, merger, by = 'name')

plot.cortex(plot.data$npos, regions = wta$label, vmin = 0, vmax = max(ptc.pos),
            cm = 'viridis', name = '') +
  labs(fill='# Positive') +
  scale_fill_colormap(colormap='viridis', limits=c(0, max(ptc.pos)), oob = scales::squish,
                      breaks = c(0, 2, 4, 6, 8, 10))
  
ggsave('oasis_ns_positivity.png', width = 8, height = 2, units = 'in')

# === Look at stage criteria met ========

ns <- pos.df
colnames(ns) <- str_replace(colnames(ns), '.WScore', '')
ns <- select(ns, all_of(ptc.order))

criteria <- data.frame(
  rid = df[df$PTCStage == 'NS', 'Subject'],
  stage1 = ns[, 1],
  stage2 = apply(ns[, 2:4], 1, any),
  stage3 = apply(ns[, 5:6], 1, any),
  stage4 = apply(ns[, 7:8], 1, any)
)

critera.yes1 <- criteria %>%
  filter(stage1 == T)

critera.not1 <- criteria %>%
  filter(stage1 == F)