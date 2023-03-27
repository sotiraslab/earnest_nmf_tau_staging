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

# === Read data ============

df <- read.csv(PATH.ADNI.STAGED) %>%
  filter(Group == 'TrainingBaseline') %>%
  mutate(PETDate = as_datetime(ymd(DateTau)))

# === Add CDR ======

cdr.adni <- cdr %>%
  dplyr::select(RID, USERDATE, CDGLOBAL) %>%
  rename(DateCDR.Long=USERDATE, CDR.Long=CDGLOBAL) %>%
  mutate(DateCDR.Long=as_datetime(ymd(DateCDR.Long)))

cdr.merged <- left_join(df, cdr.adni, by='RID') %>%
  drop_na(CDR.Long) %>%
  filter(as.Date(DateCDR.Long) >= as.Date(CDRDate)) %>%
  group_by(RID) %>%
  mutate(TickerCDR = as.integer(difftime(DateCDR.Long, first(CDRDate), units='days'))) %>%
  filter(first(CDR.Long) < 1) %>%
  ungroup() %>%
  mutate(CDR1.0Event = ifelse(CDR.Long >= 1, 1, 0))

# print subjects
unique.subs <- group_by(cdr.merged, RID) %>%
  slice_head(n = 1)

table(unique.subs$PTCStage)

# === Survival: CDR 1.0, omit NS =======

colors = c('0' = '#0072B2', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00', 'NS' = 'black')

survfit2(Surv(TickerCDR, CDR1.0Event) ~ PTCStage, data=cdr.merged) %>% 
  ggsurvfit() +
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

fit <- surv_fit(Surv(TickerCDR, CDR1.0Event) ~ PTCStage, data=cdr.merged)
surv_pvalue(fit)
pairwise_survdiff(Surv(TickerCDR, CDR1.0Event) ~ PTCStage, data=cdr.merged)
ggsave('adni_survival_cdr1.png', width=8, height=6)
