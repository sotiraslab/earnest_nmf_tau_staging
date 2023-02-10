# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.ADNI.WSCORES <- '../fig3/adni_data_with_wscores.csv'
PATH.ADNI.ORDER <- '../fig3/adni_wscore_stage_order.csv'
PATH.SCRIPT.BARPLOT <- '../../scripts/stacked_barplot.R'
PATH.SCRIPT.STAGING <- '../../scripts/stage_assigner.R'


# === Apply staging to all ADNI data ========

df <- read.csv(PATH.ADNI.WSCORES)

# Create stage column for each region
source(PATH.SCRIPT.STAGING)

stage.order <- read.csv(PATH.ADNI.ORDER)$Region

bin.w <- df %>%
  dplyr::select(contains('.WScore')) %>%
  mutate(across(.fns=function(x) ifelse(x >= 2.5, 1, 0)))

colnames(bin.w) <- gsub('.WScore', '', colnames(bin.w))
bin.w <- bin.w[, stage.order]

# see bootstrapping of W score data for stage definition
df$PTCStage <- assign.stages(bin.w, stage.order, c(1, 2, 2, 2, 3, 3, 4, 4), p='any', atypical = 'NS')

df.all <- df %>%
  arrange(RID, DateTau) %>%
  group_by(RID) %>%
  mutate(VisitNumber = rank(DateTau)) %>%
  ungroup()

# === prepare for plotting ======

# select only training data
df <- df.all %>%
  filter(Group == 'TrainingBaseline')

# define colors
stage.colors <- c('0' = 'white', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00',
                  'NS' = 'gray')

# load plotting function
source(PATH.SCRIPT.BARPLOT)

# === CDR ========

stacked.barplot(df, 'CDRBinned', 'PTCStage', colors=stage.colors) + xlab('CDR')
ggsave('adni_cdr_bar.png', width=4, height=8)

# ==== Binned Centiloid ========

df$CentiloidBinned <- cut(df$Centiloid, c(-Inf, 40, 60, 80, 100, Inf),
                          labels=c('<40', '40-60', '60-80', '80-100', '>100'))

stacked.barplot(df, 'CentiloidBinned', 'PTCStage', colors=stage.colors) + xlab('Centiloid')
ggsave('adni_centiloid_bar.png', width=6, height=8) 

# === APOE =========

stacked.barplot(df, 'HasE4', 'PTCStage', levels = c(T, F), colors=stage.colors) +
  xlab('APOE') +
  scale_x_discrete(labels=c('E4-', 'E4+'))

ggsave('adni_apoe_bar.png', width=4, height=8)

# ======= save =====

write.csv(df.all, 'adni_data_with_staging.csv', row.names = F, quote=F, na="")
