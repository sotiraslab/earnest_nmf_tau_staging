# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ggsignif))
sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.OASIS.WSCORES <- '../../derivatives/oasis3/data_with_wscores.csv'
PATH.ADNI.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'
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
df <- df.all %>%
  filter(Group == 'TrainingSet')

# define colors
stage.colors <- c('0' = 'white', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00',
                  'NS' = 'gray')

# load plotting function
source(PATH.SCRIPT.BARPLOT)

# === stage by CDR status =======

df$CDR <- factor(df$CDR, labels=c('0.0', '0.5', '1.0+'))

cdr.colors = c('0.0'='white', '0.5'='#0072B2', '1.0+'='#CC79A7')

stacked.barplot(df, 'PTCStage', 'CDR', colors=cdr.colors) +
  xlab('Stage')
ggsave('oasis_cdr_bar.png', width=6, height=8) 

data <- stacked.barplot(df, 'PTCStage', 'CDR', return.data = T)
write.csv(data, 'oasis_cdr_bar.csv')

# # ==== Bin Centiloid ========

df$CentiloidBinned <- cut(df$Centiloid, c(-Inf, 40, 60, 80, 100, Inf),
                                labels=c('<40', '40-60', '60-80', '80-100', '>100'))

stacked.barplot(df, 'CentiloidBinned', 'PTCStage', colors=stage.colors) +
  xlab('Centiloid') + 
  guides(fill = guide_legend(title='Stage'))

ggsave('oasis_centiloid_bar.png', width=6, height=8)

data <- stacked.barplot(df, 'CentiloidBinned', 'PTCStage', return.data = T)
write.csv(data, 'oasis_centiloid_bar.csv')

# === plot APOE =========

stacked.barplot(df, 'HasE4', 'PTCStage', levels=c(T, F), colors=stage.colors) +
  xlab('APOE') +
  guides(fill = guide_legend(title='Stage')) + 
  scale_x_discrete(labels=c('E4-', 'E4+'))

ggsave('oasis_apoe_bar.png', width=4, height=8)

data <- stacked.barplot(df, 'HasE4', 'PTCStage', return.data = T)
write.csv(data, 'oasis_apoe_bar.csv')

# === MMSE ===========

# run stats, get results and convert to arguments understood by ggsignif
anova <- aov(PACC.Original ~ PTCStage, data = df)
posthoc <- as.data.frame(TukeyHSD(anova, method='fdr')$PTCStage)

posthoc.res <- posthoc %>%
  rownames_to_column('comparison') %>%
  mutate(annotation = cut(`p adj`,
                          breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c('***', "**", "*", ""),
                          include.lowest = T)
  )
posthoc.sig <- filter(posthoc.res, `p adj` < 0.05)

comparisons <- str_split(posthoc.sig$comparison, '-')
n.sig <- nrow(posthoc.sig)

mean.pacc <- group_by(df, PTCStage) %>%
  summarise(PACC.Original = mean(PACC.Original, na.rm=T)) %>%
  mutate(x=seq(0.5, by=1, length.out=6),
         xend=seq(1.5, by=1, length.out=6))

ggplot(data = df, aes(x = PTCStage, y = PACC.Original, fill = PTCStage)) +
  geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
  geom_segment(data=mean.pacc, aes(x=x, xend=xend, y=PACC.Original, yend=PACC.Original),
               color='black',
               linewidth=1) + 
  scale_fill_manual(values=stage.colors) +
  theme_light() +
  xlab("Stage") +
  theme(legend.position = 'none',
        text = element_text(size=20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  coord_cartesian(ylim = c(-10, 6)) +
  scale_y_continuous(breaks=c(0, -2.5, -5, -7.5, -10)) +
  geom_signif(comparisons=comparisons,
              annotations = posthoc.sig$annotation,
              y_position = c(2, 3, 4, 5, 1, 2),
              tip_length = 0.01,
              size=.75,
              textsize = 7) +
  ylab('PACC')

ggsave('oasis_pacc_scatter.png', width=6, height=8)

print(summary(anova))

# save posthoc stats
posthoc.res <- posthoc.res %>%
  mutate(across(where(is.numeric), round, 3))
write.csv(posthoc.res, 'SUPPLEMENT_pacc_posthoc_oasis3.csv')

# save ========

path.out <- '../../derivatives/oasis3/data_with_staging.csv'
write.csv(df.all, path.out, row.names = F, quote=F, na="")
