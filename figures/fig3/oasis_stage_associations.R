# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.OASIS.WSCORES <- '../fig2/oasis_data_with_wscores.csv'
PATH.ADNI.ORDER <- '../fig2/adni_wscore_stage_order.csv'
PATH.SCRIPT.BARPLOT <- '../../scripts/stacked_barplot.R'
PATH.SCRIPT.STAGING <- '../../scripts/stage_assigner.R'

# === read data =========

# adrc data with ADNI Wscores
df.all <- read.csv(PATH.OASIS.WSCORES)

# === Create stage column for each region ===========

source(PATH.SCRIPT.STAGING)

stage.order <- read.csv(PATH.ADNI.ORDER)$Region

bin.w <- df.all %>%
  dplyr::select(contains('.WScore')) %>%
  mutate(across(.fns=function(x) ifelse(x >= 2.5, 1, 0)))

colnames(bin.w) <- gsub('.WScore', '', colnames(bin.w))
bin.w <- bin.w[, stage.order]

# see bootstrapping of W score data for stage definition
df.all$PTCStage <- assign.stages(bin.w, stage.order, c(1, 2, 2, 2, 3, 3, 4, 4), p='any', atypical = 'NS')

# === prepare for plotting ======

# select only training data
df.train <- df.all %>%
  filter(Group == 'TrainingSet')

# define colors
stage.colors <- c('0' = 'white', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00',
                  'NS' = 'gray')

# load plotting function
source(PATH.SCRIPT.BARPLOT)

# === stage by CDR status =======

df.train$CDR <- factor(df.train$CDR, labels=c('0.0', '0.5', '1.0+'))

cdr.colors = c('0.0'='white', '0.5'='#0072B2', '1.0+'='#CC79A7')

stacked.barplot(df.train, 'PTCStage', 'CDR', colors=cdr.colors) +
  xlab('Stage')
ggsave('oasis_cdr_bar.png', width=6, height=8) 

# # ==== Bin Centiloid ========

df.train$CentiloidBinned <- cut(df.train$Centiloid, c(-Inf, 40, 60, 80, 100, Inf),
                                labels=c('<40', '40-60', '60-80', '80-100', '>100'))

stacked.barplot(df.train, 'CentiloidBinned', 'PTCStage', colors=stage.colors) + xlab('Centiloid')
ggsave('oasis_centiloid_bar.png', width=6, height=8)

# === plot APOE =========

stacked.barplot(df.train, 'HasE4', 'PTCStage', levels=c(T, F), colors=stage.colors) +
  xlab('APOE') +
  scale_x_discrete(labels=c('E4-', 'E4+'))

ggsave('oasis_apoe_bar.png', width=4, height=8)

# === MMSE ===========

ggplot(data = df.train, aes(x = PTCStage, y = MMSE, fill = PTCStage)) + 
  geom_boxplot() +
  scale_fill_manual(values=stage.colors) +
  theme_light() +
  xlab("Stage") +
  theme(legend.position = 'none',
        text = element_text(size=20)) 

ggsave('oasis_mmse_boxplot.png', width=6, height=8)


# save ========
write.csv(df.all, 'oasis_data_with_staging.csv', row.names = F, quote=F, na="")
