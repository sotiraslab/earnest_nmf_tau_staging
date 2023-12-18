# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(colormap))
sh(library(ppcor))
sh(library(R.matlab))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Required Files ============

PATH.DATA.ROIS <- '../../derivatives/adni/main_data.csv'
PATH.DATA.STAGING <- '../../derivatives/adni/data_with_staging.csv'
PATH.NMF.REGIONS <- '../../derivatives/adni/nmf_regions.csv'
PATH.GGSEG <- '../../scripts/ggseg_plots.R'

# === Read data ====

df <- read.csv(PATH.DATA.ROIS) %>%
  filter(Group == 'TrainingBaseline')

tmp <- read.csv(PATH.DATA.STAGING) %>%
  filter(Group == 'TrainingBaseline')

df$PTCStage <- tmp$PTCStage

regions <- read.csv(PATH.NMF.REGIONS)$Feature

# === add amyloid rois ========

av45 <- ucberkeleyav45 %>%
  select(RID, EXAMDATE, all_of(regions)) %>%
  rename(DateAmyloid=EXAMDATE) %>%
  mutate(DateAmyloid = as.character(DateAmyloid))

fbb <- ucberkeleyfbb %>%
  select(RID, EXAMDATE, all_of(regions)) %>%
  rename(DateAmyloid=EXAMDATE) %>%
  mutate(DateAmyloid = as.character(DateAmyloid))

abeta.av45 <- df %>%
  filter(AmyloidTracer == 'AV45') %>%
  select(RID, DateAmyloid) %>%
  left_join(av45, by = c('RID', 'DateAmyloid'))

abeta.fbb <- df %>%
  filter(AmyloidTracer == 'FBB') %>%
  select(RID, DateAmyloid) %>%
  left_join(fbb, by = c('RID', 'DateAmyloid'))

abeta.big <- df %>%
  select(RID, AmyloidTracer, PTCStage) %>%
  left_join(rbind(abeta.av45, abeta.fbb), by = c('RID')) %>%
  arrange(RID)

# === create statistics ===========

# parameters
base.stage <- '0'
compare.stages <- c('1', '2', '3', '4', 'NS')
logp.cap <- 10
n.comparisons <- length(compare.stages) * length(regions)

# tau
tau.statistics <- data.frame(stage=rep(NA, n.comparisons),
                             region=rep(NA, n.comparisons),
                             tval=rep(NA, n.comparisons),
                             pval=rep(NA, n.comparisons))

for (i in seq_along(compare.stages)) {
  for (j in seq_along(regions)) {
    stage <- compare.stages[i]
    region <- regions[j]
    row <- (i - 1) * length(regions) + j
    x <- df[df$PTCStage == '0', region]
    y <- df[df$PTCStage == stage, region]
    test <- t.test(x, y, alternative = 'l')
    tau.statistics[row, ] <- c(stage, region, test$statistic, test$p.value)
  }
}


tau.statistics <- tau.statistics %>%
  mutate(pval.adjust = p.adjust(pval, method = 'fdr'),
         logp = -log10(pval.adjust),
         logp = ifelse(logp >= logp.cap, logp.cap, logp),
         logp = ifelse(pval.adjust >= 0.05, NA, logp))

# Amyloid
amy.statistics <- data.frame(stage=rep(NA, n.comparisons),
                             region=rep(NA, n.comparisons),
                             tval=rep(NA, n.comparisons),
                             pval=rep(NA, n.comparisons))

for (i in seq_along(compare.stages)) {
  for (j in seq_along(regions)) {
    stage <- compare.stages[i]
    region <- regions[j]
    row <- (i - 1) * length(regions) + j
    x <- abeta.big[abeta.big$PTCStage == '0', region]
    y <- abeta.big[abeta.big$PTCStage == stage, region]
    test <- t.test(x, y, alternative = 'l')
    amy.statistics[row, ] <- c(stage, region, test$statistic, test$p.value)
  }
}


amy.statistics <- amy.statistics %>%
  mutate(pval.adjust = p.adjust(pval, method = 'fdr'),
         logp = -log10(pval.adjust),
         logp = ifelse(logp >= logp.cap, logp.cap, logp),
         logp = ifelse(pval.adjust >= 0.05, NA, logp))

# GM volume
gm.regions <- str_replace(regions, '_SUVR', '_VOLUME')
gm.statistics <- data.frame(stage=rep(NA, n.comparisons),
                            region=rep(NA, n.comparisons),
                            tval=rep(NA, n.comparisons),
                            pval=rep(NA, n.comparisons))

for (i in seq_along(compare.stages)) {
  for (j in seq_along(gm.regions)) {
    stage <- compare.stages[i]
    region <- gm.regions[j]
    row <- (i - 1) * length(gm.regions) + j
    x <- df[df$PTCStage == '0', region]
    y <- df[df$PTCStage == stage, region]
    test <- t.test(x, y, alternative = 'g')
    gm.statistics[row, ] <- c(stage, region, test$statistic, test$p.value)
  }
}


gm.statistics <- gm.statistics %>%
  mutate(pval.adjust = p.adjust(pval, method = 'fdr'),
         logp = -log10(pval.adjust),
         logp = ifelse(logp >= logp.cap, logp.cap, logp),
         logp = ifelse(pval.adjust >= 0.05, NA, logp))

# ==== source plot script =======

source(PATH.GGSEG)
ggseg.regions <- adni.labels.to.ggseg(regions)

# ==== plot: tau =================


for (x in compare.stages) {
  plot.data <- tau.statistics %>%
    filter(stage == x)
  
  plot.cortex(values = plot.data$logp,
              regions = ggseg.regions, legend = T,
              vmin = 0, 
              vmax = logp.cap,
              cm = 'plasma', 
              name = sprintf('Stage %s Tau', x)) +
    labs(fill="-log10(p)")
  
  ggsave(sprintf('adni_tau_stage_%s.png', x), width = 8, height = 2, units = 'in')
}

# ==== plot: amyloid =================


for (x in compare.stages) {
  plot.data <- amy.statistics %>%
    filter(stage == x)
  
  plot.cortex(values = plot.data$logp,
              regions = ggseg.regions, legend = T,
              vmin = 0, 
              vmax = logp.cap,
              cm = 'bathymetry', 
              name = sprintf('Stage %s Amyloid', x)) +
    labs(fill="-log10(p)")
  
  ggsave(sprintf('adni_amyloid_stage_%s.png', x), width = 8, height = 2, units = 'in')
}

# ==== plot: GM =================

for (x in compare.stages) {
  plot.data <- gm.statistics %>%
    filter(stage == x)
  
  plot.cortex(values = plot.data$logp,
              regions = ggseg.regions, legend = T,
              vmin = 0, 
              vmax = logp.cap,
              cm = 'magma', 
              name = sprintf('Stage %s Atrophy', x)) +
    labs(fill="-log10(p)")
  
  ggsave(sprintf('adni_atrophy_stage_%s.png', x), width = 8, height = 2, units = 'in')
}

