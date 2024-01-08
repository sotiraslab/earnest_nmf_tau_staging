# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
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

# PTC Staging
stage.order <- read.csv(PATH.ADNI.ORDER)$Region
bin.w <- df.all %>%
  dplyr::select(contains('.WScore')) %>%
  mutate(across(.fns=function(x) ifelse(x >= 2.5, 1, 0)))

colnames(bin.w) <- gsub('.WScore', '', colnames(bin.w))
bin.w <- bin.w[, stage.order]
df.all$PTCStage <- assign.stages(bin.w, stage.order, c(1, 2, 2, 2, 3, 3, 4, 4), p='any', atypical = 'NS')

# Braak staging
bin.w <- df.all %>%
  dplyr::select(contains('.W') & contains("BRAAK")) %>%
  mutate(across(.fns=function(x) ifelse(x >= 2.5, 1, 0)))
colnames(bin.w) <- gsub('.W', '', colnames(bin.w))
df.all <- df.all %>%
  mutate(BraakStage = assign.stages(bin.w, colnames(bin.w), c(1, 2, 3, 4, 5), p='any', atypical = 'NS'),
         BraakStage = recode(BraakStage, '1'='I', '2'='III', '3'='IV', '4'='V', '5'='VI'),
         BraakStage = factor(BraakStage, levels=c('0', 'I', 'III', 'IV', 'V', 'VI', 'NS')))

# === prepare for plotting ======

# select only training data
df <- df.all %>%
  filter(Group == 'TrainingSet')

# define colors
stage.colors <- c('0' = 'white', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00',
                  'NS' = 'gray')
braak.colors <- c('0' = 'white', 'I'= '#009E73', 'III' = '#F0E442', 'IV' = '#E69F00', 'V' = '#D55E00',
                  'VI' = '#CC79A7', 'NS' = 'gray')


# load plotting function
source(PATH.SCRIPT.BARPLOT)

# === compare PTC & Braak assignments ========

tbl.data <- df %>%
  drop_na(BraakStage)

tbl.data <- as.data.frame(table(tbl.data$PTCStage, tbl.data$BraakStage))
colnames(tbl.data) <- c('PTCStage', 'BraakStage', 'Count')
tbl.data$BraakStage <- factor(tbl.data$BraakStage,
                              levels = c('NS', 'VI', 'V', 'IV', 'III', 'I', '0'))
tbl.data$Over <- ifelse(tbl.data$Count >= 50, "Y", 'N')

ggplot() +
  geom_tile(data=tbl.data, aes(x=PTCStage, y=BraakStage, fill=Count)) + 
  scale_fill_colormap(colormap='freesurface-blue', reverse = T,
                      limits=c(0, 50),
                      oob=scales::squish) + 
  coord_equal() +
  scale_x_discrete(expand=expansion(mult = c(0, 0))) +
  scale_y_discrete(expand=expansion(mult = c(0, 0))) + 
  geom_text(data=tbl.data, aes(x=PTCStage, y=BraakStage, label=Count, color=Over)) +
  scale_color_manual(values=c('Y'='white', 'N'='black')) +
  guides(color='none') + 
  theme(text = element_text(size=15)) + 
  ylab('Braak Stage') +
  xlab('PTC Stage')

ggsave('SUPPLEMENT_oasis_ptc_v_braak_staging.png', width = 6, height = 6, units = 'in')


# === stage by CDR status =======

df$CDR <- factor(df$CDR, labels=c('0.0', '0.5', '1.0+'))
cdr.colors = c('0.0'='white', '0.5'='#0072B2', '1.0+'='#CC79A7')

# PTC staging
stacked.barplot(df, 'PTCStage', 'CDR', colors=cdr.colors) +
  xlab('Stage')
ggsave('oasis_cdr_bar.png', width=6, height=8) 
data <- stacked.barplot(df, 'PTCStage', 'CDR', return.data = T)
write.csv(data, 'oasis_cdr_bar.csv')

# Braak staging
stacked.barplot(df, 'BraakStage', 'CDR', colors=cdr.colors, dropna=T) +
  xlab('Stage')
ggsave('BRAAK_oasis_cdr_bar.png', width=6, height=8) 
data <- stacked.barplot(df, 'BraakStage', 'CDR', return.data = T, dropna=T)
write.csv(data, 'BRAAK_oasis_cdr_bar.csv')

# # ==== Bin Centiloid ========

df$CentiloidBinned <- cut(df$Centiloid, c(-Inf, 40, 60, 80, 100, Inf),
                                labels=c('<40', '40-60', '60-80', '80-100', '>100'))

# PTC staging
stacked.barplot(df, 'CentiloidBinned', 'PTCStage', colors=stage.colors) +
  xlab('Centiloid') + 
  guides(fill = guide_legend(title='Stage'))
ggsave('oasis_centiloid_bar.png', width=6, height=8)
data <- stacked.barplot(df, 'CentiloidBinned', 'PTCStage', return.data = T)
write.csv(data, 'oasis_centiloid_bar.csv')

# Braak staging
stacked.barplot(df, 'CentiloidBinned', 'BraakStage', colors=braak.colors, dropna = T) +
  xlab('Centiloid') + 
  guides(fill = guide_legend(title='Stage'))
ggsave('BRAAK_oasis_centiloid_bar.png', width=6, height=8)
data <- stacked.barplot(df, 'CentiloidBinned', 'BraakStage', return.data = T, dropna = T)
write.csv(data, 'BRAAK_oasis_centiloid_bar.csv')

# === plot APOE =========

# PTC staging
stacked.barplot(df, 'HasE4', 'PTCStage', levels=c(T, F), colors=stage.colors) +
  xlab('APOE') +
  guides(fill = guide_legend(title='Stage')) + 
  scale_x_discrete(labels=c('E4-', 'E4+'))
ggsave('oasis_apoe_bar.png', width=4, height=8)
data <- stacked.barplot(df, 'HasE4', 'PTCStage', return.data = T)
write.csv(data, 'oasis_apoe_bar.csv')

# Braak staging
stacked.barplot(df, 'HasE4', 'BraakStage', levels=c(T, F), colors=braak.colors, dropna = T) +
  xlab('APOE') +
  guides(fill = guide_legend(title='Stage')) + 
  scale_x_discrete(labels=c('E4-', 'E4+'))
ggsave('BRAAK_oasis_apoe_bar.png', width=4, height=8)
data <- stacked.barplot(df, 'HasE4', 'BraakStage', return.data = T, dropna = T)
write.csv(data, 'BRAAK_oasis_apoe_bar.csv')

# === Chi squared =======

cramers.v <- function(chi) {
  data <- chi$observed
  k <- min(nrow(data), ncol(data)) - 1
  x <- unname(chi$statistic)
  n <- sum(data)
  v <- sqrt((x/n) / k)
  return (v)
}

# PTC Staging
chis <- list(
  'cdr' = chisq.test(table(df$PTCStage, df$CDR)),
  'centiloid' = chisq.test(table(df$PTCStage, df$CentiloidBinned)),
  'apoe' = chisq.test(table(df$PTCStage, df$HasE4))
)
chi.df <- sapply(chis, function(x) {
  c('chi'=unname(x$statistic),
    'df'=unname(x$parameter),
    'p'=round(unname(x$p.value), 3))
})
chi.df <- as.data.frame(t(chi.df))
chi.df$cramerv <- sapply(chis, cramers.v)
write.csv(chi.df, 'chi_squared_results_oasis3.csv')

# Braak Staging
chis <- list(
  'cdr' = chisq.test(table(df$BraakStage, df$CDR)),
  'centiloid' = chisq.test(table(df$BraakStage, df$CentiloidBinned)),
  'apoe' = chisq.test(table(df$BraakStage, df$HasE4))
)
chi.df <- sapply(chis, function(x) {
  c('chi'=unname(x$statistic),
    'df'=unname(x$parameter),
    'p'=round(unname(x$p.value), 3))
})
chi.df <- as.data.frame(t(chi.df))
chi.df$cramerv <- sapply(chis, cramers.v)
write.csv(chi.df, 'BRAAK_chi_squared_results_oasis3.csv')

# === PACC - PTC staging ===========

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
print(etaSquared(anova))

# save posthoc stats
posthoc.res <- posthoc.res %>%
  mutate(across(where(is.numeric), round, 3))
write.csv(posthoc.res, 'SUPPLEMENT_pacc_posthoc_oasis3.csv')

# === PACC - BRAAK staging ===========

df.braak <- df %>% drop_na(BraakStage)

# run stats, get results and convert to arguments understood by ggsignif
anova <- aov(PACC.Original ~ BraakStage, data = df.braak)
posthoc <- as.data.frame(TukeyHSD(anova, method='fdr')$BraakStage)

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

mean.pacc <- group_by(df.braak, BraakStage) %>%
  summarise(PACC.Original = mean(PACC.Original, na.rm=T)) %>%
  mutate(x=seq(0.5, by=1, length.out=7),
         xend=seq(1.5, by=1, length.out=7))

ggplot(data = df.braak, aes(x = BraakStage, y = PACC.Original, fill = BraakStage)) +
  geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
  geom_segment(data=mean.pacc, aes(x=x, xend=xend, y=PACC.Original, yend=PACC.Original),
               color='black',
               linewidth=1) + 
  scale_fill_manual(values=braak.colors) +
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

ggsave('BRAAK_oasis_pacc_scatter.png', width=6, height=8)

print(summary(anova))
print(etaSquared(anova))

# save posthoc stats
posthoc.res <- posthoc.res %>%
  mutate(across(where(is.numeric), round, 3))
write.csv(posthoc.res, 'BRAAK_pacc_posthoc_oasis3.csv')

# ===== Stage proportions by Centiloid ========

cent.data <- df %>%
  select(Centiloid, PTCStage) %>%
  arrange(Centiloid) %>%
  mutate(TmpStage = as.numeric(ifelse(PTCStage == 'NS', 0, PTCStage)),
         Stage1 = cumsum(TmpStage >= 1) / n(),
         Stage2 = cumsum(TmpStage >= 2) / n(),
         Stage3 = cumsum(TmpStage >= 3) / n(),
         Stage4 = cumsum(TmpStage >= 4) / n(),
         StageNS = cumsum(PTCStage == 'NS') / n()) %>%
  pivot_longer(c(Stage1, Stage2, Stage3, Stage4, StageNS), names_to = 'Stage', values_to = 'Proportion') %>%
  mutate(Stage = gsub("Stage", '', Stage))

ggplot(cent.data, aes(x=Centiloid, y=Proportion, color=Stage, linetype=Stage)) +
  geom_step(linewidth=1) +
  theme_light() +
  scale_color_manual(values=unname(stage.colors[c('1', '2', '3', '4', 'NS')])) +
  coord_cartesian(ylim=c(0, .5)) +
  theme(text = element_text(size=20))

ggsave('SUPPLEMENT_oasis_stage_by_centiloids.png', width=8, height=6)

# save ========

path.out <- '../../derivatives/oasis3/data_with_staging.csv'
write.csv(df.all, path.out, row.names = F, quote=F, na="")
