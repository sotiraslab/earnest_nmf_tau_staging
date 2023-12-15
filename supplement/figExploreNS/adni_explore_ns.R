# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.DATA <- '../../derivatives/adni/data_with_staging.csv'
PATH.PTC.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'

# === Read ========

df <- read.csv(PATH.DATA)
ptc.order <- read.csv(PATH.PTC.ORDER)$Region

df <- df %>%
  filter(Group == 'TrainingBaseline') %>%
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
    fml <- as.formula(sprintf('%s ~ Stageable + CorticalTauAverage', dependent))
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
                  'PACC.ADNI',
                  'MMSE',
                  'Centiloid',
                  'CorticalTauAverage',
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
write.csv(stat.data, 'stats_raw.csv', quote = F, na = '', row.names = F)

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

write.csv(str.table, 'stats_formatted.csv', quote = F, na = '', row.names = F)

# === Old functions ========

# t.test.statisticss <- function(dependent) {
#   mask.stageable <- df$Stageable == 'Stageable'
#   mask.ns <- df$Stageable == 'NS'
#   mask.0 <- df$Stageable == 'Stage 0'
#   
#   data.stageable <- df[mask.stageable, dependent]
#   data.ns <- df[mask.ns, dependent]
#   data.0 <- df[mask.0, dependent]
#   
#   # 0 vs stageable 
#   t1 <- t.test(x = data.0,
#                y = data.stageable,
#                alternative = 't')
#   
#   
#   # 0 vs NS
#   t2 <- t.test(x = data.0,
#                y = data.ns,
#                alternative = 't')
#   
#   # NS vs stageable
#   t3 <- t.test(x = data.ns,
#                y = data.stageable,
#                alternative = 't')
#   
#   
#   vals <- c(dependent,
#             mean(data.0, na.rm = T),
#             sd(data.0, na.rm = T),
#             mean(data.stageable, na.rm = T),
#             sd(data.stageable, na.rm = T),
#             mean(data.ns, na.rm = T),
#             sd(data.ns, na.rm = T),
#             t1$p.value,
#             t2$p.value,
#             t3$p.value)
#   
# }
# 
# lm.statistics <- function(dependent, adjust.tau = F) {
#   mask.stageable <- df$Stageable == 'Stageable'
#   mask.ns <- df$Stageable == 'NS'
#   mask.0 <- df$Stageable == 'Stage 0'
#   
#   # modeling
#   data1 <- df[mask.0 | mask.stageable, ]
#   data2 <- df[mask.0 | mask.ns, ]
#   data3 <- df[mask.ns | mask.stageable, ]
#   
#   if (adjust.tau) {
#     fml <- as.formula(sprintf('%s ~ Stageable + CorticalTauAverage', dependent))
#   } else {
#     fml <- as.formula(sprintf('%s ~ Stageable', dependent))
#   }
#   
#   m1 <- lm(fml, data = data1)
#   m2 <- lm(fml, data = data2)
#   m3 <- lm(fml, data = data3)
#   
#   data.stageable <- df[mask.stageable, dependent]
#   data.ns <- df[mask.ns, dependent]
#   data.0 <- df[mask.0, dependent]
#   
#   vals <- c(dependent,
#             mean(data.0, na.rm = T),
#             sd(data.0, na.rm = T),
#             mean(data.stageable, na.rm = T),
#             sd(data.stageable, na.rm = T),
#             mean(data.ns, na.rm = T),
#             sd(data.ns, na.rm = T),
#             summary(m1)$coefficients[2, 4],
#             summary(m2)$coefficients[2, 4],
#             summary(m3)$coefficients[2, 4])
# }
# 
# 
# get.stats <- function() {
#   dependents <- c('Age',
#                   'PACC.ADNI',
#                   'MMSE',
#                   'Centiloid',
#                   ptc.order,
#                   'CorticalTauAverage',
#                   'Laterality')
#   
#   mat <- matrix(NA, nrow = length(dependents), ncol = 10)
#   t.data <- as.data.frame(mat)
#   colnames(t.data) <- c('Variable',
#                         'Stage0.mean',
#                         'Stage0.sd',
#                         'Stageable.mean',
#                         'Stageable.sd',
#                         'NS.mean',
#                         'NS.sd',
#                         'Stage0.vs.Stageable',
#                         'Stage0.vs.NS',
#                         'NS.vs.Stageable')
#   
#   for (i in seq_along(dependents)) {
#     dependent <- dependents[i]
#     t.data[i, ] <- get.stats.for.dependent(dependent)
#   }
#   return (t.data)
# }

