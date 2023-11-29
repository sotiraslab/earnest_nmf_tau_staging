# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
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

# PTC staging
# see bootstrapping of W score data for stage definition
bin.w <- df %>%
  dplyr::select(contains('.WScore')) %>%
  mutate(across(.fns=function(x) ifelse(x >= 2.5, 1, 0)))
colnames(bin.w) <- gsub('.WScore', '', colnames(bin.w))
bin.w <- bin.w[, stage.order]
df$PTCStage <- assign.stages(bin.w, stage.order, c(1, 2, 2, 2, 3, 3, 4, 4), p='any', atypical = 'NS')

# Braak staging
bin.w <- df %>%
  dplyr::select(contains('.W') & contains("BRAAK")) %>%
  mutate(across(.fns=function(x) ifelse(x >= 2.5, 1, 0)))
colnames(bin.w) <- gsub('.W', '', colnames(bin.w))
df <- df %>%
  mutate(BraakStage = assign.stages(bin.w, colnames(bin.w), c(1, 2, 3, 4, 5), p='any', atypical = 'NS'),
         BraakStage = recode(BraakStage, '1'='I', '2'='III', '3'='IV', '4'='V', '5'='VI'),
         BraakStage = factor(BraakStage, levels=c('0', 'I', 'III', 'IV', 'V', 'VI', 'NS')))

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
tbl.data$Over100 <- ifelse(tbl.data$Count >= 100, "Y", 'N')

ggplot() +
  geom_tile(data=tbl.data, aes(x=PTCStage, y=BraakStage, fill=Count)) + 
  scale_fill_colormap(colormap='freesurface-blue', reverse = T,
                      limits=c(0, 100),
                      oob=scales::squish) + 
  coord_equal() +
  scale_x_discrete(expand=expansion(mult = c(0, 0))) +
  scale_y_discrete(expand=expansion(mult = c(0, 0))) + 
  geom_text(data=tbl.data, aes(x=PTCStage, y=BraakStage, label=Count, color=Over100)) +
  scale_color_manual(values=c('Y'='white', 'N'='black')) +
  guides(color='none') + 
  theme(text = element_text(size=15)) + 
  ylab('Braak Stage') +
  xlab('PTC Stage')

ggsave('SUPPLEMENT_adni_ptc_v_braak_staging.png', width = 6, height = 6, units = 'in')

# === CDR ========

cdr.colors = c('0.0'='white', '0.5'='#0072B2', '1.0+'='#CC79A7')

# ptc staging
stacked.barplot(df, 'PTCStage', 'CDRBinned', colors=cdr.colors) +
  xlab('Stage') +
  guides(fill=guide_legend(title='CDR'))
ggsave('adni_cdr_bar.png', width=6, height=8)
data <- stacked.barplot(df, 'PTCStage', 'CDRBinned', return.data = T)
write.csv(data, 'adni_cdr_bar.csv')

# braak staging
stacked.barplot(df, 'BraakStage', 'CDRBinned', colors=cdr.colors) +
  xlab('Stage') +
  guides(fill=guide_legend(title='CDR'))
ggsave('BRAAK_adni_cdr_bar.png', width=6, height=8)
data <- stacked.barplot(df, 'BraakStage', 'CDRBinned', return.data = T)
write.csv(data, 'BRAAK_adni_cdr_bar.csv')

# ==== Binned Centiloid ========

df$CentiloidBinned <- cut(df$Centiloid, c(-Inf, 40, 60, 80, 100, Inf),
                          labels=c('<40', '40-60', '60-80', '80-100', '>100'))

# ptc staging
stacked.barplot(df, 'CentiloidBinned', 'PTCStage', colors=stage.colors) +
  xlab('Centiloid') +
  guides(fill = guide_legend(title='Stage'))
ggsave('adni_centiloid_bar.png', width=6, height=8) 
data <- stacked.barplot(df, 'PTCStage', 'CentiloidBinned', return.data = T)
write.csv(data, 'adni_centiloid_bar.csv')

# braak staging
stacked.barplot(df, 'CentiloidBinned', 'BraakStage', colors=braak.colors) +
  xlab('Centiloid') +
  guides(fill = guide_legend(title='Stage'))
ggsave('BRAAK_adni_centiloid_bar.png', width=6, height=8) 
data <- stacked.barplot(df, 'BraakStage', 'CentiloidBinned', return.data = T)
write.csv(data, 'BRAAK_adni_centiloid_bar.csv')

# === APOE =========

# ptc staging
stacked.barplot(df, 'HasE4', 'PTCStage', levels = c(T, F), colors=stage.colors) +
  xlab('APOE') +
  guides(fill = guide_legend(title='Stage')) + 
  scale_x_discrete(labels=c('E4-', 'E4+'))
ggsave('adni_apoe_bar.png', width=4, height=8)
data <- stacked.barplot(df, 'HasE4', 'PTCStage', return.data = T)
write.csv(data, 'adni_apoe_bar.csv')

# braak staging
stacked.barplot(df, 'HasE4', 'BraakStage', levels = c(T, F), colors=braak.colors) +
  xlab('APOE') +
  guides(fill = guide_legend(title='Stage')) + 
  scale_x_discrete(labels=c('E4-', 'E4+'))
ggsave('BRAAK_adni_apoe_bar.png', width=4, height=8)
data <- stacked.barplot(df, 'HasE4', 'BraakStage', return.data = T)
write.csv(data, 'BRAAK_adni_apoe_bar.csv')

# === Chi squared =======

# ptc staging
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
write.csv(chi.df, 'chi_squared_results_adni.csv')

# braak staging
chis <- list(
  'cdr' = chisq.test(table(df$BraakStage, df$CDRBinned)),
  'centiloid' = chisq.test(table(df$BraakStage, df$CentiloidBinned)),
  'apoe' = chisq.test(table(df$BraakStage, df$HasE4))
)
chi.df <- sapply(chis, function(x) {
  c('chi'=unname(x$statistic),
    'df'=unname(x$parameter),
    'p'=round(unname(x$p.value), 3))
})
chi.df <- as.data.frame(t(chi.df))
write.csv(chi.df, 'BRAAK_chi_squared_results_adni.csv')

# === PACC - PTC Staging =========

# run stats, get results and convert to arguments understood by ggsignif
anova <- aov(PACC.ADNI ~ PTCStage, data = df)
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
  summarise(PACC.ADNI = mean(PACC.ADNI, na.rm=T)) %>%
  mutate(x=seq(0.5, by=1, length.out=6),
         xend=seq(1.5, by=1, length.out=6))

ggplot(data = df, aes(x = PTCStage, y = PACC.ADNI, fill = PTCStage)) +
  geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
  geom_segment(data=mean.pacc, aes(x=x, xend=xend, y=PACC.ADNI, yend=PACC.ADNI),
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
              y_position = c(1, 2, 3, 4, 5, 1, 2),
              tip_length = 0.01,
              size=.75,
              textsize = 7) +
  ylab('PACC')

ggsave('adni_pacc_scatter.png', width=6, height=8)

print(summary(anova))

# save posthoc stats
posthoc.res <- posthoc.res %>%
  mutate(across(where(is.numeric), round, 3))
write.csv(posthoc.res, 'SUPPLEMENT_pacc_posthoc_adni.csv')

# === PACC - BRAAK Staging =========

anova <- aov(PACC.ADNI ~ BraakStage, data = df)
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

mean.pacc <- group_by(df, BraakStage) %>%
  summarise(PACC.ADNI = mean(PACC.ADNI, na.rm=T)) %>%
  mutate(x=seq(0.5, by=1, length.out=7),
         xend=seq(1.5, by=1, length.out=7))

ggplot(data = df, aes(x = BraakStage, y = PACC.ADNI, fill = BraakStage)) +
  geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
  geom_segment(data=mean.pacc, aes(x=x, xend=xend, y=PACC.ADNI, yend=PACC.ADNI),
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
              y_position = c(1, 2, 3, 4, 5, 1, 2),
              tip_length = 0.01,
              size=.75,
              textsize = 7) +
  ylab('PACC')

ggsave('BRAAK_adni_pacc_scatter.png', width=6, height=8)

print(summary(anova))

# save posthoc stats
posthoc.res <- posthoc.res %>%
  mutate(across(where(is.numeric), round, 3))
write.csv(posthoc.res, 'BRAAK_pacc_posthoc_adni.csv')

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

ggsave('SUPPLEMENT_adni_stage_by_centiloids.png', width=8, height=6)

# ======= save =====

path.out <- '../../derivatives/adni/data_with_staging.csv'
write.csv(df.all, path.out, row.names = F, quote=F, na="")
