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

PATH.OASIS.STAGED <- '../fig4/oasis_data_with_staging.csv'
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

ggsave('oasis_survival_cdr1.png', width=8, height=6)
pairwise_survdiff(Surv(Ticker, CDR1.0Event) ~ PTCStage, data=cdr.merged)

# === Survival: CDR 1.0, omit NS =======

cdr.merged.no.NS <- filter(cdr.merged, PTCStage != 'NS')

# eight stage group
survfit2(Surv(Ticker, CDR1.0Event) ~ PTCStage, data=cdr.merged.no.NS) %>% 
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

ggsave('oasis_survival_cdr1_omit_NS.png', width=8, height=6)
pairwise_survdiff(Surv(Ticker, CDR1.0Event) ~ PTCStage, data=cdr.merged.no.NS)
