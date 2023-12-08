
# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Files needed ===========

PATH.FTP <- '../../rawdata/oasis_flortaucipir.csv'
PATH.CENTILOID <- '../../derivatives/oasis3/manual_centiloid.csv'
PATH.CLINICAL <- '../../rawdata/OASIS3_data_files/scans/UDSb4-Form_B4__Global_Staging__CDR__Standard_and_Supplemental/resources/csv/files/OASIS3_UDSb4_cdr.csv'
PATH.DEMO <- '../../rawdata/OASIS3_data_files/scans/demo-demographics/resources/csv/files/OASIS3_demographics.csv'
PATH.NEUROPSYCH <- '../../rawdata/OASIS3_data_files/scans/pychometrics-Form_C1__Cognitive_Assessments/resources/csv/files/OASIS3_UDSc1_cognitive_assessments.csv'
PATH.OASIS.SUBS <- '../../subject_ids/oasis_subjects.csv'
PATH.PACC.SCRIPT <- '../../scripts/pacc.R'
PATH.FS <- '../../rawdata/oasis_freesurfer.csv'

# === Load tau ROI data =============

# This contains all the ROI SUVR values extracted from
# PUP.  Both PVC and non-PVC.  This at least gives
# the very maximum number of subjects that can be used
# for NMF (434).  But need to filter out cases that are
# amyloid negative.

big.df <- read.csv(PATH.FTP)

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

amyloid.df <- read.csv(PATH.CENTILOID)


amyloid.df.proc <- amyloid.df %>%
  select(PUP_PUPTIMECOURSEDATA.ID,
         tracer,
         ManualCentiloid) %>%
  rename(AmyloidTracer = tracer,
         Centiloid = ManualCentiloid,
         PUPLongID = PUP_PUPTIMECOURSEDATA.ID) %>%
  mutate(Subject = str_extract(PUPLongID, 'OAS\\d+'),
         SessionAmyloid = adrc_session_to_number(PUPLongID),
         AmyloidPositive = ifelse(AmyloidTracer == 'AV45',
                                  Centiloid >= 20.6,
                                  Centiloid >= 16.4),
         Indexer = 1:n())

amyloid.selector <- function(Indexer, TauAmyloidGap, AmyloidStatus) {
  
  omit1 <- is.na(TauAmyloidGap) 
  omit2 <- (TauAmyloidGap < -365)
  mask.omit <- omit1 | omit2 
  
  if (all(mask.omit)) {
    return(NA)
  }
  
  Indexer <- Indexer[ ! mask.omit]
  TauAmyloidGap <- TauAmyloidGap[! mask.omit]
  AmyloidStatus <- AmyloidStatus[! mask.omit]
  
  closest.scan.index <- which.min(abs(TauAmyloidGap))
  
  if (abs(TauAmyloidGap[closest.scan.index]) <= 365 | AmyloidStatus[closest.scan.index] == 1) {
    return (Indexer[which.min(abs(TauAmyloidGap))])
  } else {
    return (NA)
  }
}

df.merge <- left_join(df, amyloid.df.proc, by='Subject') %>%
  mutate(TauAmyloidDiff = Session - SessionAmyloid) %>%
  group_by(Subject) %>%
  summarise(Indexer=amyloid.selector(Indexer, TauAmyloidDiff, AmyloidPositive)) %>%
  ungroup() %>%
  left_join(select(amyloid.df.proc, -Subject), by='Indexer')

# assign status
df <- left_join(df, df.merge, by='Subject') %>%
  mutate(TauAmyloidDiff = Session - SessionAmyloid)

# === Add CDR ==========

clinical.df <- read.csv(PATH.CLINICAL)

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
  mutate(CDR = ifelse(TauClinicalDiff > 365, NA, CDR)) %>%
  ungroup()

# === Add Demographics ==========

demo <- read.csv(PATH.DEMO)

demo.proc <- demo %>%
  select(OASISID, AgeatEntry, AgeatDeath, GENDER, APOE) %>%
  rename(Subject=OASISID) %>%
  mutate(HasE4 = grepl('4', APOE))

df <- left_join(df, demo.proc, by='Subject')
df$Age <- df$AgeatEntry + (df$Session / 365.25)

# === Add MMSE ==========

mmse <- read.csv(PATH.CLINICAL)

mmse <- mmse %>%
  select(OASISID, MMSE, days_to_visit) %>%
  rename(SessionMMSE=days_to_visit,
         Subject=OASISID) %>%
  select(Subject, MMSE, SessionMMSE)

df <- left_join(df, mmse, by='Subject') %>%
  mutate(TauMMSEDiff = abs(Session - SessionMMSE)) %>%
  group_by(Subject) %>%
  slice_min(TauMMSEDiff, with_ties = F) %>%
  ungroup()

# === Add neuropsych ==========

nps <- read.csv(PATH.NEUROPSYCH)

nps.cols <- c('LOGIMEM',
              'MEMUNITS',
              'digsym',
              'srtfree',
              'asscmem',
              'ANIMALS',
              'VEG',
              'tma',
              'tmb',
              'bnt')

nps <- nps %>%
  select(OASISID, days_to_visit, all_of(nps.cols)) %>%
  rename(SessionNeuroPsych = days_to_visit,
         Subject=OASISID)

df <- left_join(df, nps, by='Subject') %>%
  mutate(TauNPSDiff = abs(Session - SessionNeuroPsych)) %>%
  group_by(Subject) %>%
  slice_min(TauNPSDiff, with_ties = F) %>%
  ungroup()

bad.nps <- (df$TauNPSDiff > 365) | (is.na(df$TauNPSDiff))
df[bad.nps, nps.cols] <- NA

# === Add regional tau =========

tau.rois <- read.csv('../../rawdata/oasis_flortaucipir.csv')

# === Select amyloid positive df ==========

df.amyloidpos <- df[! is.na(df$AmyloidPositive) & df$AmyloidPositive, ]
table(df.amyloidpos$CDR, useNA = 'always')

# corpus callosum & cerebllum are omitted to match ADNI
pat.inc <- "(PET_fSUVR_[LR]_CTX_.*)"
pat.exc <- "CRPCLM|CBLL"
cols <- colnames(tau.rois)[grepl(pat.inc, colnames(tau.rois), perl = T) &
                           ! grepl(pat.exc, colnames(tau.rois), perl = T)]

# ROIS from ADRC have different names than ADNI/ggseg
# but are ordered the same
adni.rois <- read.csv('../../derivatives/adni/nmf_regions.csv')
converter <-  data.frame(ADRC=cols, ADNI=adni.rois$Feature)

tau.rois$Subject <- str_extract(tau.rois$PUP_PUPTIMECOURSEDATA.ID, 'OAS\\d+')

# add total cortical mean for a course tau index
joiner <- tau.rois[, c("Subject", cols, 'PET_fSUVR_TOT_CORTMEAN',
                       'PET_fSUVR_L_AMYGDALA',
                       'PET_fSUVR_L_HIPPOCAMPUS',
                       'PET_fSUVR_R_AMYGDALA',
                       'PET_fSUVR_R_HIPPOCAMPUS')]
colnames(joiner) <- c("Subject", converter$ADNI, 'TotalCtxTauMean',
                      'LEFT_AMYGDALA_SUVR',
                      'LEFT_HIPPOCAMPUS_SUVR',
                      'RIGHT_AMYGDALA_SUVR',
                      'RIGHT_HIPPOCAMPUS_SUVR')

df <- left_join(df, joiner, by='Subject')

# === add Volumes =======

# read volumes
volume.df <-  read.csv(PATH.FS)
all.cols <- colnames(volume.df)
vol.cols <- all.cols[str_detect(all.cols, '^(lh|rh).*_volume$')]
vol.cols <- c(vol.cols,
              'Left.Amygdala_volume',
              'Left.Hippocampus_volume',
              'Right.Amygdala_volume',
              'Right.Hippocampus_volume')

vol.edit <- volume.df[, c('FS_FSDATA.ID', vol.cols)]
colnames(vol.edit) <- c('FSId', toupper(vol.cols))

df <- left_join(df, vol.edit, by='FSId')

# === add Braak regions ======

volume.weighted.mean <- function(pet.data, vol.data, search.columns) {
  
  pattern <- paste(search.columns, collapse='|')
  pet.cols <- colnames(pet.data)[str_detect(colnames(pet.data), pattern)]
  vol.cols <- colnames(vol.data)[str_detect(colnames(vol.data), pattern)]
  
  pet.data <- pet.data[, pet.cols]
  vol.data <- vol.data[, vol.cols]
  print(dim(pet.data))
  print(dim(vol.data))
  
  volumes.norm <- vol.data / rowSums(vol.data)
  pet.norm <- pet.data * volumes.norm
  result <- rowSums(pet.norm)
  
  return(result)
}

braak1.regs <- c('ENTORHINAL')

braak3.regs <- c('PARAHIPPOCAMPAL',
                 'FUSIFORM',
                 'LINGUAL',
                 'AMYGDALA')

braak4.regs <- c('MIDDLETEMPORAL',
                 'CAUDALANTERIORCINGULATE',
                 'ROSTRALANTERIORCINGULATE',
                 'POSTERIORCINGULATE',
                 'ISTHMUSCINGULATE',
                 'INSULA',
                 'INFERIORTEMPORAL',
                 'TEMPORALPOLE')

braak5.regs <- c('SUPERIORFRONTAL',
                 'LATERALORBITOFRONTAL',
                 'MEDIALORBITOFRONTAL',
                 'FRONTALPOLE',
                 'CAUDALMIDDLEFRONTAL',
                 'ROSTRALMIDDLEFRONTAL',
                 'PARSOPERCULARIS',
                 'PARSORBITALIS',
                 'PARSTRIANGULARIS',
                 'LATERALOCCIPITAL',
                 'SUPRAMARGINAL',
                 'INFERIORPARIETAL',
                 'SUPERIORTEMPORAL',
                 'SUPERIORPARIETAL',
                 'PRECUNEUS',
                 'BANKSSTS',
                 'TRANSVERSETEMPORAL')

braak6.regs <- c('PERICALCARINE',
                 'POSTCENTRAL',
                 'CUNEUS',
                 'PRECENTRAL',
                 'PARACENTRAL')

df.pet <- df %>%
  select(matches('_SUVR'))
df.vol <- df %>%
  select(matches('_VOLUME'))

df$BRAAK1_SUVR <- volume.weighted.mean(df.pet, df.vol, braak1.regs)
df$BRAAK3_SUVR <- volume.weighted.mean(df.pet, df.vol, braak3.regs)
df$BRAAK4_SUVR <- volume.weighted.mean(df.pet, df.vol, braak4.regs)
df$BRAAK5_SUVR <- volume.weighted.mean(df.pet, df.vol, braak5.regs)
df$BRAAK6_SUVR <- volume.weighted.mean(df.pet, df.vol, braak6.regs)

# === Add groups =============

df$Group <- NA

df$Group <- ifelse(is.na(df$Group) & df$Subject %in% df.amyloidpos$Subject,
                   'TrainingSet',
                   NA)

cdr0 <- (df$CDR == 0 & ! is.na(df$CDR))
amyneg <- (! is.na(df$AmyloidPositive) & df$AmyloidPositive == F)
both <- cdr0 & amyneg
df$Group <- ifelse(is.na(df$Group) & both,
                   'ControlSet',
                   df$Group)

df$Group <- ifelse(is.na(df$Group),
                   'Other',
                   df$Group)

# === Calculate PACC ======

source(PATH.PACC.SCRIPT)

df$PACC.Original <- compute.pacc(df,
                                 pacc.columns = c('srtfree', 'MEMUNITS', 'digsym', 'MMSE'),
                                 cn.mask <- df$Group == 'ControlSet',
                                 higher.better = c(T, T, T, T))

# === Compute NPS composites =========

# not the pacc, but can still use "compute.pacc" function

df$Composite.MEM <- compute.pacc(df,
                                 pacc.columns = c('LOGIMEM', 'srtfree', 'asscmem'),
                                 cn.mask = df$Group == 'ControlSet',
                                 higher.better = c(T, T, T),
                                 min.required = 2)

df$Composite.EF <- compute.pacc(df,
                                pacc.columns = c('tma', 'tmb'),
                                cn.mask = df$Group == 'ControlSet',
                                higher.better = c(F, F),
                                min.required = 2)

df$Composite.LANG <- compute.pacc(df,
                                  pacc.columns = c('ANIMALS', 'VEG', 'bnt'),
                                  cn.mask = df$Group == 'ControlSet',
                                  higher.better = c(T, T, T, T),
                                  min.required = 2)

# === Generate subject list =========

# this code creates the subject-ids file for OASIS3
# NOTE: overwriting the subject lists in the subject_ids folder
#       may result in a different set of subjects being analyzed!

# sub.ids <- select(df, Subject, Session, Group)
# write.csv(sub.ids, 'oasis_subjects.csv', row.names = F)


# === Apply Subject list =========

original.subs <- read.csv(PATH.OASIS.SUBS)
current.subs <- select(df, -Group)
df.original.subs <- left_join(original.subs, current.subs, by=c('Subject', 'Session'))

# === Save ==========

write.csv(df.original.subs, '../../derivatives/oasis3/main_data.csv', na='', quote=F, row.names = F)
  