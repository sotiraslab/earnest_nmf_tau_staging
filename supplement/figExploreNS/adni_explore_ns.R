# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.DATA <- '../../derivatives/adni/data_with_staging.csv'

# === Read ========

df <- read.csv(PATH.DATA)

df <- df %>%
  filter(Group == 'TrainingBaseline') %>%
  mutate(Stageable = ifelse(PTCStage %in% c('1', '2', '3', '4'), 'Stageable', NA),
         Stageable = ifelse(PTCStage == '0', 'Stage 0', Stageable),
         Stageable = ifelse(PTCStage == 'NS', 'NS', Stageable))

# === Create stats DF ========

get.stats.for.dependent <- function(dependent) {
  mask.stageable <- df$Stageable == 'Stageable'
  mask.ns <- df$Stageable == 'NS'
  mask.0 <- df$Stageable == 'Stage 0'
  
  data.stageable <- df[mask.stageable, dependent]
  data.ns <- df[mask.ns, dependent]
  data.0 <- df[mask.0, dependent]
  
  # 0 vs stageable 
  t1 <- t.test(x = data.0,
               y = data.stageable,
               alternative = 't')
  
  
  # 0 vs NS
  t2 <- t.test(x = data.0,
               y = data.ns,
               alternative = 't')
  
  # NS vs stageable
  t3 <- t.test(x = data.ns,
               y = data.stageable,
               alternative = 't')
  
  
  vals <- c(dependent,
            mean(data.0, na.rm = T),
            sd(data.0, na.rm = T),
            mean(data.stageable, na.rm = T),
            sd(data.stageable, na.rm = T),
            mean(data.ns, na.rm = T),
            sd(data.ns, na.rm = T),
            t1$p.value,
            t2$p.value,
            t3$p.value)
  
}

get.stats <- function() {
  dependents <- c('Age',
                  'PACC.ADNI',
                  'MMSE',
                  'Centiloid',
                  'CorticalTauAverage')
  
  mat <- matrix(NA, nrow = length(dependents), ncol = 10)
  t.data <- as.data.frame(mat)
  colnames(t.data) <- c('Variable',
                        'Stage0.mean',
                        'Stage0.sd',
                        'Stageable.mean',
                        'Stageable.sd',
                        'NS.mean',
                        'NS.sd',
                        'Stage0.vs.Stageable',
                        'Stage0.vs.NS',
                        'NS.vs.Stageable')
  
  for (i in seq_along(dependents)) {
    dependent <- dependents[i]
    t.data[i, ] <- get.stats.for.dependent(dependent)
  }
  return (t.data)
}

t.data <- get.stats()


