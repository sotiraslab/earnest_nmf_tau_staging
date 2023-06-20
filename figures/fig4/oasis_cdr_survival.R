# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ggsurvfit))
sh(library(lubridate))
sh(library(survival))
sh(library(survminer))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Necessary Files =========

PATH.OASIS.STAGED <- '../../derivatives/oasis3/data_with_staging.csv'
PATH.OASIS.CDR <- '../../rawdata/OASIS3_data_files/scans/UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv'

# === Read data =======

df <- read.csv(PATH.OASIS.STAGED) %>%
  filter(Group == 'TrainingSet')

# === Add longitudinal CDR ======

oasis_session_to_number <- function(col) {
  extracted <- str_extract(col, 'd\\d+')
  no.leading.d <- substr(extracted, 2, nchar(extracted))
  number <- as.numeric(no.leading.d)
  
  return (number)
}

cdr.oasis <- read.csv(PATH.OASIS.CDR) %>%
  select(OASIS_session_label, CDRTOT) %>%
  mutate(Subject=str_extract(OASIS_session_label, 'OAS\\d+'),
         SessionCDR.Long=oasis_session_to_number(OASIS_session_label)) %>%
  rename(CDR.Long=CDRTOT) %>%
  select(-OASIS_session_label)

cdr.merged <- left_join(df, cdr.oasis, by='Subject') %>%
  drop_na(CDR.Long) %>%
  filter(SessionCDR.Long >= SessionCDR) %>%
  group_by(Subject) %>%
  filter(first(CDR.Long) < 1) %>%
  mutate(Ticker = SessionCDR.Long - SessionCDR) %>%
  ungroup() %>%
  mutate(CDR1.0Event = ifelse(CDR.Long >= 1, 1, 0),
         CDR0.5Event = ifelse(CDR.Long >= 0.5, 1, 0))

# see supplement for plot with NS included
cdr.merged <- filter(cdr.merged, PTCStage != 'NS')

# print subjects
unique.subs <- group_by(cdr.merged, Subject) %>%
  slice_head(n = 1)

table(unique.subs$PTCStage)
  
# === Survival: CDR 1.0, with NS =======

colors = c('0' = '#0072B2', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00', 'NS' = 'gray')

# eight stage group
survfit2(Surv(Ticker, CDR1.0Event) ~ PTCStage, data=cdr.merged) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Probability of CDR<1"
  ) +
  add_confidence_interval() +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  theme_ggsurvfit_KMunicate() +
  theme(text = element_text(size=20),
        legend.position = 'bottom')

fit <- surv_fit(Surv(Ticker, CDR1.0Event) ~ PTCStage, data=cdr.merged)
posthoc <- pairwise_survdiff(Surv(Ticker, CDR1.0Event) ~ PTCStage, data=cdr.merged)
ggsave('oasis_survival_cdr1.png', width=8, height=6)

surv_pvalue(fit)

# ==========

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
write.csv(ps, 'SUPPLEMENT_survival_posthoc_oasis3.csv')
