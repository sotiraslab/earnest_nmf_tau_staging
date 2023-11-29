# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(ggsurvfit))
sh(library(lubridate))
sh(library(survival))
sh(library(survminer))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Necessary Files =========

PATH.ADNI.STAGED <- '../../derivatives/adni/data_with_staging.csv'
PATH.EXAMDATE <- '../../scripts/adni_examdate.R'

# === Read data ============

df <- read.csv(PATH.ADNI.STAGED) %>%
  filter(Group == 'TrainingBaseline') %>%
  mutate(PETDate = as_datetime(ymd(DateTau)))

# === Add CDR ======

source(PATH.EXAMDATE)

cdr.adni <- cdr %>%
  mutate(DateCDR.Long = ifelse(is.na(EXAMDATE),
                               get.examdate.from.registry(cdr),
                               EXAMDATE),
         DateCDR.Long = as_datetime(ymd(DateCDR.Long))) %>%
  dplyr::select(RID, DateCDR.Long, CDGLOBAL) %>%
  rename(CDR.Long=CDGLOBAL)

cdr.merged <- left_join(df, cdr.adni, by='RID') %>%
  drop_na(CDR.Long) %>%
  filter(as.Date(DateCDR.Long) >= as.Date(CDRDate)) %>%
  group_by(RID) %>%
  mutate(TickerCDR = as.integer(difftime(DateCDR.Long, first(CDRDate), units='days'))) %>%
  filter(first(CDR.Long) < 1) %>%
  ungroup() %>%
  mutate(CDR1.0Event = ifelse(CDR.Long >= 1, 1, 0))

# see supplement for plot with NS included
cdr.merged.ptc <- filter(cdr.merged, PTCStage != 'NS')
cdr.merged.braak <- filter(cdr.merged, BraakStage != 'NS')

# print subjects
unique.subs <- group_by(cdr.merged, RID) %>%
  slice_head(n = 1)

table(unique.subs$PTCStage)

# === Survival: CDR 1.0, omit NS - PTC Staging =======

colors = c('0' = '#0072B2', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00', 'NS' = 'gray')

survfit2(Surv(TickerCDR, CDR1.0Event) ~ PTCStage, data=cdr.merged.ptc) %>% 
  ggsurvfit(linetype_aes = T, size=1.5) +
  labs(
    x = "Days",
    y = "Probability of CDR<1",
  ) +
  add_confidence_interval() +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  theme_ggsurvfit_KMunicate() +
  theme(text = element_text(size=20),
        legend.position = 'bottom')

# stats, including NS
fit <- surv_fit(Surv(TickerCDR, CDR1.0Event) ~ PTCStage, data=cdr.merged.ptc)
posthoc <- pairwise_survdiff(Surv(TickerCDR, CDR1.0Event) ~ PTCStage, data=cdr.merged.ptc)
ggsave('adni_survival_cdr1.png', width=8, height=6)

surv_pvalue(fit)

ps <- posthoc$p.value %>%
  as.data.frame() %>%
  rownames_to_column('a') %>%
  pivot_longer(cols = c('0', '1', '2', '3'), values_to = 'p.value', names_to = 'b') %>%
  mutate(comparison = paste(a, b, sep=' vs ')) %>%
  select(comparison, p.value) %>%
  filter(! is.na(p.value)) %>%
  mutate(annotation = cut(p.value,
                          breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c('***', "**", "*", ""),
                          include.lowest = T),
         p.value = round(p.value, 5))
write.csv(ps, 'SUPPLEMENT_survival_posthoc_adni.csv', row.names = F)

# === Survival: CDR 1.0, omit NS - Braak Staging =======

colors = c('0' = '#0072B2', 'I'= '#009E73', 'III' = '#F0E442', 'IV' = '#E69F00', 'V' = '#D55E00',
           'VI' = '#CC79A7', 'NS' = 'gray')

survfit2(Surv(TickerCDR, CDR1.0Event) ~ BraakStage, data=cdr.merged.braak) %>% 
  ggsurvfit(linetype_aes = T, size=1.5) +
  labs(
    x = "Days",
    y = "Probability of CDR<1",
  ) +
  add_confidence_interval() +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  theme_ggsurvfit_KMunicate() +
  theme(text = element_text(size=20),
        legend.position = 'bottom')

# stats, including NS
fit <- surv_fit(Surv(TickerCDR, CDR1.0Event) ~ BraakStage, data=cdr.merged.braak)
posthoc <- pairwise_survdiff(Surv(TickerCDR, CDR1.0Event) ~ BraakStage, data=cdr.merged.braak)
ggsave('BRAAK_adni_survival_cdr1.png', width=8, height=6)

surv_pvalue(fit)

ps <- posthoc$p.value %>%
  as.data.frame() %>%
  rownames_to_column('a') %>%
  pivot_longer(cols = c('0', 'I', 'III', 'IV', 'V'), values_to = 'p.value', names_to = 'b') %>%
  mutate(comparison = paste(a, b, sep=' vs ')) %>%
  select(comparison, p.value) %>%
  filter(! is.na(p.value)) %>%
  mutate(annotation = cut(p.value,
                          breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c('***', "**", "*", ""),
                          include.lowest = T),
         p.value = round(p.value, 5))
write.csv(ps, 'BRAAK_survival_posthoc_adni.csv', row.names = F)
