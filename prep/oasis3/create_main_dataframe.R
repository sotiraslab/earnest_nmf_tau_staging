
# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Load tau ROI data =============

# This contains all the ROI SUVR values extracted from
# PUP.  Both PVC and non-PVC.  This at least gives
# the very maximum number of subjects that can be used
# for NMF (434).  But need to filter out cases that are
# amyloid negative.

path <- '../../rawdata/oasis_flortaucipir.csv'
big.df <- read.csv(path)

# extract only the subjects
adrc_session_to_number <- function(col) {
  extracted <- str_extract(col, 'd\\d+')
  no.leading.d <- substr(extracted, 2, nchar(extracted))
  number <- as.numeric(no.leading.d)
  
  return (number)
}

df <- big.df[, c("PUP_PUPTIMECOURSEDATA.ID", 'FSId', 'MRId'), drop=0]
df$Subject <- str_extract(df$PUP_PUPTIMECOURSEDATA.ID, 'OAS\\d+')
df$Session <- adrc_session_to_number(df$PUP_PUPTIMECOURSEDATA.ID)

# === Add Amyloid Status ===========

# based on the amyloid cortical values
# see check_centiloid_conversion.R for a little more explanation

path <- '../../derivatives/oasis3/manual_centiloid.csv'
amyloid.df <- read.csv(path)


amyloid.df.proc <- amyloid.df %>%
  select(PUP_PUPTIMECOURSEDATA.ID,
         tracer,
         ManualCentiloid) %>%
  rename(AmyloidTracer = tracer,
         CentiloidPVC = ManualCentiloid,
         PUPLongID = PUP_PUPTIMECOURSEDATA.ID) %>%
  mutate(Subject = str_extract(PUPLongID, 'OAS\\d+'),
         SessionAmyloid = adrc_session_to_number(PUPLongID),
         AmyloidPositive = ifelse(AmyloidTracer == 'AV45',
                                  CentiloidPVC >= 20.6,
                                  CentiloidPVC >= 16.4))

# merge into the df

total.tau <- nrow(df)
tau.no.amyloid <- length(df$Subject[! df$Subject %in% amyloid.df.proc$Subject])

df.merge <- left_join(df, amyloid.df.proc, by='Subject') %>%
  mutate(TauAmyloidDiff = Session - SessionAmyloid) %>%
  drop_na(AmyloidPositive)

n1 <- length(unique(df.merge$Subject))
df.merge <- filter(df.merge, TauAmyloidDiff > -365)
n2 <- length(unique(df.merge$Subject))
no.recent.amyloid <- n1 - n2
  
status <- df.merge %>%
  group_by(Subject) %>%
  arrange(-CentiloidPVC) %>%
  summarise(AnyAmyloidPositive = any(AmyloidPositive),
            MixStatus = any(AmyloidPositive) & ! all(AmyloidPositive),
            Centiloid = first(CentiloidPVC),
            TauAmyloidDiff = first(TauAmyloidDiff)) %>%
  ungroup()

# report the availability
cat("Tau subjects/images:", total.tau, '\n')
cat("With any amyloid imaging:", total.tau - tau.no.amyloid, '\n')
cat("With contemporaneous amyloid imaging:", total.tau - tau.no.amyloid - no.recent.amyloid, '\n')
cat("Amyloid+ subjects:", sum(status$AnyAmyloidPositive), '\n')
cat("Subjects with mixed amyloid status:", sum(status$MixStatus), '\n')

# assign status
df <- left_join(df, status, by='Subject')

# === Add CDR ==========

path <- '../../rawdata/OASIS3_data_files/scans/UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv'
clinical.df <- read.csv(path)

clinical.df.proc <- clinical.df %>%
  select(OASIS_session_label, CDRTOT, CDRSUM) %>%
  mutate(Subject=str_extract(OASIS_session_label, 'OAS\\d+'),
         SessionCDR=adrc_session_to_number(OASIS_session_label)) %>%
  rename(CDR=CDRTOT, CDRSumBoxes=CDRSUM) %>%
  select(-OASIS_session_label)

df <- left_join(df, clinical.df.proc, by='Subject') %>%
  mutate(TauClinicalDiff = abs(Session - SessionCDR)) %>%
  group_by(Subject) %>%
  slice_min(TauClinicalDiff, with_ties = F) %>%
  mutate(CDR = ifelse(TauClinicalDiff > 365, NA, CDR))

# === Add Demographics ==========

path <- '../../rawdata/OASIS3_data_files/scans/demo-demographics/resources/csv/files/OASIS3_demographics.csv'
demo <- read.csv(path)

demo.proc <- demo %>%
  select(OASISID, AgeatEntry, AgeatDeath, GENDER, APOE) %>%
  rename(Subject=OASISID) %>%
  mutate(HasE4 = grepl('4', APOE))

df <- left_join(df, demo.proc, by='Subject')
df$Age <- df$AgeatEntry + (df$Session / 365.25)

# === Add MMSE ==========

path <- '../../rawdata/OASIS3_data_files/scans/UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv'
mmse <- read.csv(path)

mmse <- mmse %>%
  select(OASISID, MMSE, days_to_visit) %>%
  rename(SessionMMSE=days_to_visit,
         Subject=OASISID) %>%
  select(Subject, MMSE, SessionMMSE)

df <- left_join(df, mmse, by='Subject') %>%
  mutate(TauMMSEDiff = abs(Session - SessionMMSE)) %>%
  group_by(Subject) %>%
  slice_min(TauMMSEDiff, with_ties = F)

# === Add regional tau =========

tau.rois <- read.csv('../../rawdata/oasis_flortaucipir.csv')

# === Select amyloid positive df ==========

df.amyloidpos <- df[! is.na(df$AnyAmyloidPositive) & df$AnyAmyloidPositive, ]
table(df.amyloidpos$CDR, useNA = 'always')

# corpus callosum & cerebllum are omitted to match ADNI
pat.inc <- "PET_fSUVR_[LR]_CTX_.*"
pat.exc <- "CRPCLM|CBLL"
cols <- colnames(tau.rois)[grepl(pat.inc, colnames(tau.rois), perl = T) &
                           ! grepl(pat.exc, colnames(tau.rois), perl = T)]

# ROIS from ADRC have different names than ADNI/ggseg
# but are ordered the same
adni.rois <- read.csv('../../derivatives/adni/nmf_regions.csv')
converter <-  data.frame(ADRC=cols, ADNI=adni.rois$Feature)

tau.rois$Subject <- str_extract(tau.rois$PUP_PUPTIMECOURSEDATA.ID, 'OAS\\d+')

# add total cortical mean for a course tau index
joiner <- tau.rois[, c("Subject", cols, 'PET_fSUVR_TOT_CORTMEAN')]
colnames(joiner) <- c("Subject", converter$ADNI, 'TotalCtxTauMean')

df <- left_join(df, joiner, by='Subject')

# === Add groups =============

df$Group <- NA

df$Group <- ifelse(is.na(df$Group) & df$Subject %in% df.amyloidpos$Subject,
                   'TrainingSet',
                   NA)

cdr0 <- (df$CDR == 0 & ! is.na(df$CDR))
amyneg <- (! is.na(df$AnyAmyloidPositive) & df$AnyAmyloidPositive == F)
both <- cdr0 & amyneg
df$Group <- ifelse(is.na(df$Group) & both,
                   'ControlSet',
                   df$Group)

df$Group <- ifelse(is.na(df$Group),
                   'Other',
                   df$Group)

# === Save ==========

write.csv(df, '../../derivatives/oasis3/main_data.csv', na='', quote=F, row.names = F)
  