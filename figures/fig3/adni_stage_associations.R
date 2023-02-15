# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ggsignif))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.ADNI.WSCORES <- '../../derivatives/adni/data_with_wscores.csv'
PATH.ADNI.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'
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

cdr.colors = c('0.0'='white', '0.5'='#0072B2', '1.0+'='#CC79A7')

stacked.barplot(df, 'PTCStage', 'CDRBinned', colors=cdr.colors) +
  xlab('Stage') +
  guides(fill=guide_legend(title='CDR'))

ggsave('adni_cdr_bar.png', width=6, height=8)

data <- stacked.barplot(df, 'PTCStage', 'CDRBinned', return.data = T)
write.csv(data, 'adni_cdr_bar.csv')

# ==== Binned Centiloid ========

df$CentiloidBinned <- cut(df$Centiloid, c(-Inf, 40, 60, 80, 100, Inf),
                          labels=c('<40', '40-60', '60-80', '80-100', '>100'))

stacked.barplot(df, 'CentiloidBinned', 'PTCStage', colors=stage.colors) + xlab('Centiloid')
ggsave('adni_centiloid_bar.png', width=6, height=8) 

data <- stacked.barplot(df, 'PTCStage', 'CentiloidBinned', return.data = T)
write.csv(data, 'adni_centiloid_bar.csv')

# === APOE =========

stacked.barplot(df, 'HasE4', 'PTCStage', levels = c(T, F), colors=stage.colors) +
  xlab('APOE') +
  scale_x_discrete(labels=c('E4-', 'E4+'))

ggsave('adni_apoe_bar.png', width=4, height=8)

data <- stacked.barplot(df, 'HasE4', 'PTCStage', return.data = T)
write.csv(data, 'adni_apoe_bar.csv')

# === Chi squared =======

chis <- list(
  'cdr' = chisq.test(table(df$PTCStage, df$CDRBinned)),
  'centiloid' = chisq.test(table(df$PTCStage, df$CentiloidBinned)),
  'apoe' = chisq.test(table(df$PTCStage, df$HasE4))
)

chi.df <- sapply(chis, function(x) {
  c('chi'=unname(x$statistic),
    'df'=unname(x$parameter),
    'p'=round(unname(x$p.value), 3))
})

chi.df <- as.data.frame(t(chi.df))
write.csv(chi.df, 'chi_squared_results.csv')

# === MMSE =========

# run stats, get results and convert to arguments understood by ggsignif
anova <- aov(MMSE ~ PTCStage, data = df)
posthoc <- as.data.frame(TukeyHSD(anova, method='fdr')$PTCStage)

posthoc.sig <- posthoc %>%
  filter(`p adj` < 0.05) %>%
  rownames_to_column('comparison') %>%
  mutate(annotation = cut(`p adj`,
                          breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c('***', "**", "*", ""),
                          include.lowest = T)
         )

comparisons <- str_split(posthoc.sig$comparison, '-')
n.sig <- nrow(posthoc.sig)

ggplot(data = df, aes(x = PTCStage, y = MMSE, fill = PTCStage)) +
  geom_boxplot() +
  scale_fill_manual(values=stage.colors) +
  theme_light() +
  xlab("Stage") +
  theme(legend.position = 'none',
        text = element_text(size=20)) +
  geom_signif(comparisons=comparisons,
              annotations = posthoc.sig$annotation,
              y_position = c(30.5, 32.5, 34.5, 32.5, 30.5),
              tip_length = 0.01)

ggsave('adni_mmse_boxplot.png', width=6, height=8)

# ======= save =====

path.out <- '../../derivatives/adni/data_with_staging.csv'
write.csv(df.all, path.out, row.names = F, quote=F, na="")
