# This script creates an initial dataframe of tau scans and associtated
# variables for ADNI.  

# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set wd ===========

setwd(this.dir())

# === Files needed ========

PATH.PTODEMOG.CSV <- '../../rawdata/PTDEMOG.csv'
PATH.SUBJECTS.ADNI <- '../../subject_ids/adni_subjects.csv'
PATH.PACC.SCRIPT <- '../../scripts/pacc.R'
PATH.EXAMDATE.SCRIPT <- '../../scripts/adni_examdate.R'

# === load the tau data ==========

tau.adni <- ucberkeleyav1451
ctx.cols <- colnames(tau.adni)[grepl('CTX_.*_SUVR', colnames(tau.adni), perl=T)]
ctx.cols <- ctx.cols[! grepl('UNKNOWN', ctx.cols)]
sbctx.cols <- colnames(tau.adni)[grepl('(AMYGDALA|HIPPOCAMPUS).*_SUVR', colnames(tau.adni), perl=T)]
vol.cols <- gsub('SUVR', 'VOLUME', c(ctx.cols, sbctx.cols))

df <- tau.adni %>%
  select(RID, EXAMDATE, all_of(ctx.cols), all_of(sbctx.cols), all_of(vol.cols), VISCODE,
         INFERIORCEREBELLUM_SUVR,META_TEMPORAL_SUVR) %>%
  rename(DateTau = EXAMDATE) %>%
  mutate(DateTau = as_datetime(ymd(DateTau)),
         ScanIndex = 1:n())

# === add Braak regions ============

# these are manually computed so the same approach can be applied in OASIS
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

# === normalize tau columns ===========

normalize.cols <- colnames(df)[grepl('SUVR|BRAAK|META_TEMPORAL', colnames(df))]

df <- df %>%
  mutate(across(all_of(normalize.cols), function(x) x / INFERIORCEREBELLUM_SUVR))

# === define cortical tau measure average ========
# use for splitting data

df$CorticalTauAverage <- rowMeans(df[, ctx.cols])

ctx.vol.cols <- gsub('SUVR', 'VOLUME', ctx.cols)
suv <- df[, ctx.cols]
vol <- df[, ctx.vol.cols ]
vol <- vol / rowSums(vol)
suv.weighted <- suv * vol
df$WeightedCorticalTauAverage <- rowSums(suv.weighted)

# === Add amyloid information =========

av45.record <- ucberkeleyav45[, c('RID', 'EXAMDATE', 'SUMMARYSUVR_WHOLECEREBNORM', 'SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF')]
colnames(av45.record) <- c("RID", "DateAmyloid", "CompositeSUVR", "AmyloidStatus")
av45.record$DateAmyloid <- as_datetime(ymd(av45.record$DateAmyloid))
av45.record$AmyloidTracer <- "AV45"

fbb.record <- ucberkeleyfbb[, c('RID', 'EXAMDATE', 'SUMMARYSUVR_WHOLECEREBNORM', 'SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF')]
colnames(fbb.record) <- c("RID", "DateAmyloid", "CompositeSUVR", "AmyloidStatus")
fbb.record$DateAmyloid <- as_datetime(ymd(fbb.record$DateAmyloid))
fbb.record$AmyloidTracer <- "FBB"

all.amy <- rbind(av45.record, fbb.record)
all.amy$Indexer <- 1:nrow(all.amy)

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

df.allamy <- left_join(df, all.amy, by='RID') %>%
  mutate(TauAmyloidGap=as.numeric(difftime(DateTau, DateAmyloid, units='days'))) %>%
  group_by(ScanIndex) %>%
  summarise(Indexer=amyloid.selector(Indexer, TauAmyloidGap, AmyloidStatus)) %>%
  ungroup() %>%
  left_join(all.amy, by='Indexer') %>%
  select(-RID) 

df <- left_join(df, df.allamy, by='ScanIndex') %>%
  mutate(TauAmyloidGap=as.numeric(difftime(DateTau, DateAmyloid, units='days')))

# ==== Add centiloids =======

fbb_to_centiloid <- function(suvr) {
  return ((159.08 * suvr) - 151.65)
}

av45_to_centiloid <- function(suvr) {
  return ((196.9 * suvr) - 196.03)
}

df$Centiloid <- NA
df$Centiloid <- ifelse(df$AmyloidTracer == 'AV45' & ! is.na(df$AmyloidTracer),
                       av45_to_centiloid(df$CompositeSUVR),
                       df$Centiloid)
df$Centiloid <- ifelse(df$AmyloidTracer == 'FBB' & ! is.na(df$AmyloidTracer),
                       fbb_to_centiloid(df$CompositeSUVR),
                       df$Centiloid)

# ==== add CDR ===========

source(PATH.EXAMDATE.SCRIPT)

thr <- 365 # how many days close the CDR needs to be

cdr.record <- cdr %>%
  as.data.frame() %>%
  mutate(CDRDate = ifelse(is.na(EXAMDATE),
                          get.examdate.from.registry(cdr),
                          EXAMDATE),
         CDRDate = as_datetime(ymd(CDRDate))) %>%
  select(RID, CDRDate, CDGLOBAL, CDRSB) %>%
  rename(CDRGlobal=CDGLOBAL, CDRSumBoxes=CDRSB) %>%
  drop_na(CDRGlobal)

cdr.df <- left_join(df, cdr.record, by='RID')
cdr.df$TauCDRGap <- as.numeric(abs(cdr.df$DateTau - cdr.df$CDRDate), units='days')

cdr.df <- group_by(cdr.df, ScanIndex) %>%
  slice_min(order_by=TauCDRGap, with_ties = F) %>% ungroup()

bad <- is.na(cdr.df$CDRGlobal) | (cdr.df$TauCDRGap > thr)
cdr.df[bad, c("CDRGlobal", "CDRSumBoxes")] <- NA

cdr.df <- cdr.df %>% mutate(CDRBinned=cut(CDRGlobal, breaks=c(0, .5, 1, Inf), right=F))
levels(cdr.df$CDRBinned) <- c("0.0", "0.5", "1.0+")

df <- as.data.frame(cdr.df)

# === add MMSE ==========

# add MMSE
mmse.adni <- mmse %>%
  mutate(DateMMSE = ifelse(is.na(EXAMDATE),
                           get.examdate.from.registry(mmse),
                           EXAMDATE),
         DateMMSE = as_datetime(ymd(DateMMSE))) %>%
  dplyr::select(RID, DateMMSE, MMSCORE) %>%
  rename(MMSE=MMSCORE)

mmse.merged <- left_join(df, mmse.adni, by='RID') %>%
  mutate(TauMMSEDiff = difftime(DateTau, DateMMSE, units='days')) %>%
  group_by(ScanIndex) %>%
  slice_min(abs(TauMMSEDiff), with_ties = F) %>%
  ungroup() %>%
  mutate(MMSE = ifelse(! is.na(TauMMSEDiff) & abs(TauMMSEDiff) > 365, NA, MMSE))  

df <- mmse.merged

# === Add APOE ==========

a1 <- select(apoeres, RID, APGEN1, APGEN2)
a2 <- select(apoego2, RID, APGEN1, APGEN2)
a3 <- select(apoe3, RID, APGEN1, APGEN2)

all.apoe <- do.call(rbind, list(a1, a2, a3))

all.apoe <- all.apoe %>%
  mutate(APOEGenotype=paste(
    pmin(all.apoe$APGEN1, all.apoe$APGEN2),
    pmax(all.apoe$APGEN1, all.apoe$APGEN2),
    sep='/')
  )

df <- left_join(df, all.apoe, by='RID')
df$HasE4 <- ifelse(is.na(df$APOEGenotype), NA, grepl('4', df$APOEGenotype))

# === Add demographics ==========

# age: find the age at earliest baseline scan
# calculate difference from that age

min.ages <- ptdemog %>%
  mutate(DateDemogBL = ifelse(is.na(EXAMDATE),
                              get.examdate.from.registry(ptdemog),
                              as.character(EXAMDATE)),
         DateDemogBL = as_datetime(ymd(DateDemogBL))) %>%
  select(RID, DateDemogBL, AGE, PTGENDER) %>%
  rename(AgeBL=AGE, Gender=PTGENDER) %>%
  mutate(AgeBL = as.numeric(AgeBL)) %>%
  drop_na(AgeBL) %>%
  group_by(RID) %>%
  slice_min(DateDemogBL)

df.age <- left_join(df, min.ages, by='RID')
df.age$TimeTauSinceBL <- as.numeric(difftime(df.age$DateTau, df.age$DateDemogBL, units='days')) / 365.25
df.age$Age <- as.numeric(df.age$AgeBL + df.age$TimeTauSinceBL)

# manually add some missing ages
# these are not in ADNIMERGE::ptdemog but are in the downloaded study tables

missing.age <- df.age[is.na(df.age$Age), ]
replace.rids <- missing.age$RID

demog.csv <- read.csv(PATH.PTODEMOG.CSV) %>%
  select(RID, PTDOBMM, PTDOBYY) %>%
  mutate(DOB=as.POSIXct(paste(PTDOBMM, 15, PTDOBYY, sep='/'), format='%m/%d/%Y')) %>%
  filter(RID %in% replace.rids) %>%
  select(RID, DOB)

df.age <- left_join(df.age, demog.csv, by='RID') %>%
  mutate(Age = ifelse(
    is.na(Age),
    as.numeric(difftime(DateTau, DOB, units='days') / 365.25),
    Age)
    )
  
df <- df.age

# same for some missing genders
df$Gender <- as.character(df$Gender)
gender.missing <- df[is.na(df$Gender), ]

demog.csv <- read.csv(PATH.PTODEMOG.CSV) %>%
  select(RID, PTGENDER) %>%
  filter(RID %in% gender.missing$RID) %>%
  mutate(Gender.Missing=recode(PTGENDER, `1`='Male', `2`='Female')) %>%
  select(-PTGENDER)

df.gender <- left_join(df, demog.csv, by='RID') %>%
  mutate(Gender = ifelse(is.na(Gender), Gender.Missing, Gender)) %>%
  select(-Gender.Missing)

df <- df.gender

# === Add neuropsych ===========

# for computation of PACC, need some neuropsych & ADAS cog Q4
# see https://adni.bitbucket.io/reference/pacc.html

nps.cols <- c('LIMMTOTAL',
              'LDELTOTAL',
              'AVDEL30MIN',
              'AVTOT1',
              'AVTOT2',
              'AVTOT3',
              'AVTOT4',
              'AVTOT5',
              'RAVLT.IMMEDIATE',
              'DSPANFOR',
              'DSPANBAC',
              'TRAASCOR',
              'TRABSCOR',
              'CATANIMSC',
              'CATVEGESC',
              'BNTTOTAL',
              'MINTTOTAL')

# 1. Neuropsych battery
nps <- neurobat %>%
  mutate(DateNeuropsych = ifelse(is.na(EXAMDATE),
                                get.examdate.from.registry(neurobat),
                                EXAMDATE),
         DateNeuropsych = as_datetime(ymd(DateNeuropsych)),
         BNTTOTAL = as.numeric(BNTTOTAL),
         CATVEGESC = as.numeric(CATVEGESC)) %>%
  select(RID, DateNeuropsych, all_of(nps.cols))
nps$RAVLT.IMMEDIATE.MANUAL <- rowSums(nps[, c('AVTOT1', 'AVTOT2', 'AVTOT3', 'AVTOT4', 'AVTOT5')])

df <- left_join(df, nps, by='RID') %>%
  mutate(DiffTauNPS = as.numeric(abs(difftime(DateTau, DateNeuropsych, units='days')))) %>%
  group_by(ScanIndex) %>%
  slice_min(DiffTauNPS, with_ties = F)

bad.nps <- (df$DiffTauNPS > thr) | (is.na(df$DiffTauNPS))
df[bad.nps, c(nps.cols, 'RAVLT.IMMEDIATE.MANUAL')] <- NA

# 2. ADAS
adascog <- adas %>%
  mutate(DateADAS = ifelse(is.na(EXAMDATE),
                           get.examdate.from.registry(adas),
                           EXAMDATE),
         DateADAS = as_datetime(ymd(DateADAS))) %>%
  select(RID, DateADAS, Q4SCORE) %>%
  rename(ADASQ4=Q4SCORE)

df <- left_join(df, adascog, by='RID') %>%
  mutate(DiffTauADAS = as.numeric(abs(difftime(DateTau, DateADAS, units='days')))) %>%
  group_by(ScanIndex) %>%
  slice_min(DiffTauNPS, with_ties = F)

bad.adas <- (df$DiffTauADAS > thr) | (is.na(df$DiffTauADAS))
df[bad.adas, c('ADASQ4')] <- NA

# === Create baseline training df =========

# Note that the order of operations is important for determining
# exactly who is included.  I am first selecting the amyloid
# positive individuals, then selecting the earliest scan per person.
# This should keep people who covert from Amyloid- to Amyloid+.

df.amyloidpos.bl <- df[! is.na(df$AmyloidStatus), ]
df.amyloidpos.bl <- df.amyloidpos.bl[df.amyloidpos.bl$AmyloidStatus == 1, ]

df.amyloidpos.bl <- group_by(df.amyloidpos.bl, RID) %>%
  slice_min(DateTau)

# === Create groupings ========== 

# The dataframe here contains all tau scans, not just the ADS set.
# Create a "Group" column to keep track of that

# 0. init column
df$Group <- NA

# 1. Mark images already used for NMF training
df$Group <- ifelse(is.na(df$Group) & df$ScanIndex %in% df.amyloidpos.bl$ScanIndex,
                    'TrainingBaseline',
                    df$Group)

# 2. Mark longitudinal data of subjects used for training
training.subs <- df[df$Group == 'TrainingBaseline' & ! is.na(df$Group), 'RID', drop=T]
df$Group <- ifelse(is.na(df$Group) & df$RID %in% training.subs,
                    'TrainingLongitudinal',
                    df$Group)

# sanity check that remaining subjects do not have any amyloid positivity
# apos should be 0
unassigned <- df[is.na(df$Group), ]
apos <- sum(unassigned$AmyloidStatus[! is.na(unassigned$AmyloidStatus)])

# 3. Select a control cohort - CDR0 and amyloid negative,
# who stay CDR0
mask <- (
  is.na(df$Group) &
    (! is.na(df$CDRBinned) & df$CDRBinned == '0.0') &
    (! is.na(df$AmyloidStatus) & df$AmyloidStatus == 0)
)
controls <- df[mask, ]
controls.bl <- group_by(controls, RID) %>%
  slice_min(DateTau) %>%
  ungroup()

df$Group <- ifelse(is.na(df$Group) & df$ScanIndex %in% controls.bl$ScanIndex,
                    'ControlBaseline',
                    df$Group)

# 4. Longitudinal data for controls
df$Group <- ifelse(is.na(df$Group) & df$RID %in% controls$RID,
                    'ControlLongitudinal',
                    df$Group)

# 5. Remainder should be amyloid status 0 or unknown
# and CDR>0 or unknown

remainder <- df[is.na(df$Group), ]
remainder.bl <- group_by(remainder, RID) %>%
  slice_min(DateTau) %>%
  ungroup()

df$Group <- ifelse(is.na(df$Group) & df$ScanIndex %in% remainder.bl$ScanIndex,
                    'RemainderBaseline',
                    df$Group)
df$Group <- ifelse(is.na(df$Group) & df$RID %in% remainder.bl$RID,
                   'RemainderLongitudinal',
                   df$Group)

# === Compute PACC =========

# done relative to baseline CN group
# this is using the modified formula recommended by ADNIMERGE R
# https://adni.bitbucket.io/reference/pacc.html

source(PATH.PACC.SCRIPT)

df$PACC.ADNI <- compute.pacc(df,
                             pacc.columns = c('ADASQ4', 'LDELTOTAL', 'TRABSCOR', 'MMSE'),
                             cn.mask = df$Group == 'ControlBaseline',
                             higher.better = c(F, T, F, T),
                             min.required = 2)

# === Compute NPS composites =========

# not the pacc, but can still use "compute.pacc" function

df$Composite.MEM <- compute.pacc(df,
                                 pacc.columns = c('LIMMTOTAL', 'LDELTOTAL', 'RAVLT.IMMEDIATE.MANUAL', 'AVDEL30MIN'),
                                 cn.mask = df$Group == 'ControlBaseline',
                                 higher.better = c(T, T, T, T),
                                 min.required = 2)

df$Composite.EF <- compute.pacc(df,
                                pacc.columns = c('TRAASCOR', 'TRABSCOR'),
                                cn.mask = df$Group == 'ControlBaseline',
                                higher.better = c(F, F),
                                min.required = 2)

df$Composite.LANG <- compute.pacc(df,
                                  pacc.columns = c('CATANIMSC', 'BNTTOTAL', 'MINTTOTAL'),
                                  cn.mask = df$Group == 'ControlBaseline',
                                  higher.better = c(T, T, T, T),
                                  min.required = 2)


# === Generate subject list =========

# this code creates the subject-ids file for ADNI
# NOTE: overwriting the subject lists in the subject_ids folder
#       may result in a different set of subjects being analyzed!

# sub.ids <- select(df, RID, VISCODE, Group)
# write.csv(sub.ids, 'adni_subjects.csv', row.names = F)

# === Apply subject list ==========

# this code tries to make sure the data created here
# has the same subjects used in the original paper

original.subs <- read.csv(PATH.SUBJECTS.ADNI)
current.subs <- select(df, -Group) %>%
  mutate(VISCODE=as.character(VISCODE))
df.original.subs <- left_join(original.subs, current.subs, by=c('RID', 'VISCODE'))

# === Save ===========

write.csv(df.original.subs, '../../derivatives/adni/main_data.csv', quote=F, na='', row.names=F)
