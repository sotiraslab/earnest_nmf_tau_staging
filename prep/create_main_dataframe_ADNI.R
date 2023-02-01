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

# === load the tau data ==========

tau.adni <- ucberkeleyav1451
ctx.cols <- colnames(tau.adni)[grepl('CTX_.*_SUVR', colnames(tau.adni), perl=T)]
ctx.cols <- ctx.cols[! grepl('UNKNOWN', ctx.cols)]
vol.cols <- gsub('SUVR', 'VOLUME', ctx.cols)

df <- tau.adni %>%
  select(RID, EXAMDATE, all_of(ctx.cols), all_of(vol.cols), VISCODE, INFERIORCEREBELLUM_SUVR,
         BRAAK1_SUVR, BRAAK34_SUVR, BRAAK56_SUVR, META_TEMPORAL_SUVR) %>%
  rename(DateTau = EXAMDATE) %>%
  mutate(DateTau = as_datetime(ymd(DateTau)),
         ScanIndex = 1:n())

# === normalize tau columns ===========

normalize.cols <- colnames(df)[grepl('SUVR|BRAAK|META_TEMPORAL', colnames(df))]

df <- df %>%
  mutate(across(all_of(normalize.cols), function(x) x / INFERIORCEREBELLUM_SUVR))

# === define cortical tau measure average ========
# use for splitting data

df$CorticalTauAverage <- rowMeans(df[, ctx.cols])

suv <- df[, ctx.cols]
vol <- df[, vol.cols]
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
df.allamy <- left_join(df, all.amy, by='RID')
df.allamy$TauAmyloidGap <- difftime(df.allamy$DateTau, df.allamy$DateAmyloid, units='days')
df.allamy <- drop_na(df.allamy, CompositeSUVR)

rows <- list()
c <- 0

for (i in 1:nrow(df)) {
  cat('Scan', i ,'...\n')
  row <- df[i, ]
  row$DateAmyloid <- NA
  row$CompositeSUVR <- NA
  row$AmyloidStatus <- NA
  row$TauAmyloidGap <- NA
  row$AmyloidTracer <- NA
  
  idx = df[i, 'ScanIndex']
  slice <- df.allamy[df.allamy$ScanIndex == idx, ]
  
  if (all(is.na(slice$TauAmyloidGap))) {
  } else if (min(abs(slice$TauAmyloidGap)) <= 365) {
    row <- slice[which.min(abs(slice$TauAmyloidGap)), ]
  } else {
    amy.before.tau <- slice[slice$TauAmyloidGap > 0, ]
    if (nrow(amy.before.tau) != 0) {
      most.recent <- amy.before.tau[which.min(amy.before.tau$TauAmyloidGap), ]
      if (most.recent[['AmyloidStatus']] == 1) row <- most.recent
    }
  }
  
  if (! is.null(row)) {
    rows[[i]] <- row
  }
}

C <- data.frame(do.call(rbind.data.frame, rows))
df <- C

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

thr <- 365 # how many days close the CDR needs to be
# this is the "1year" of "nmf_on_rois_final_1year"

cdr.record <- cdr
cdr.record <- cdr[, c("RID", "USERDATE", "CDGLOBAL", "CDRSB")]
colnames(cdr.record) <- c("RID", "CDRDate", "CDRGlobal", "CDRSumBoxes")
cdr.record$CDRDate <- as_datetime(ymd(cdr.record$CDRDate))

cdr.df <- left_join(df, cdr.record, by='RID')
cdr.df$TauCDRGap <- as.numeric(abs(cdr.df$DateTau - cdr.df$CDRDate), units='days')

cdr.df <- group_by(cdr.df, ScanIndex) %>%
  slice_min(order_by=TauCDRGap, with_ties = F) %>% ungroup()

bad <- is.na(cdr.df$CDRGlobal) | (cdr.df$TauCDRGap > thr)
cdr.df[bad, c("CDRDate", "CDRGlobal", "CDRSumBoxes", "TauCDRGap")] <- NA

cdr.df <- cdr.df %>% mutate(CDRBinned=cut(CDRGlobal, breaks=c(0, .5, 1, Inf), right=F))
levels(cdr.df$CDRBinned) <- c("0.0", "0.5", "1.0+")

df <- as.data.frame(cdr.df)

# === add MMSE ==========

# add MMSE
mmse.adni <- mmse %>%
  dplyr::select(RID, USERDATE, MMSCORE) %>%
  rename(DateMMSE=USERDATE, MMSE=MMSCORE) %>%
  mutate(DateMMSE=as_datetime(ymd(DateMMSE)))

mmse.merged <- left_join(df, mmse.adni, by='RID') %>%
  mutate(TauMMSEDiff = difftime(DateTau, DateMMSE, units='days')) %>%
  group_by(ScanIndex) %>%
  slice_min(abs(TauMMSEDiff), with_ties = F) %>%
  ungroup() %>%
  mutate(MMSE = ifelse(! is.na(TauMMSEDiff) & abs(TauMMSEDiff) > 365, NA, MMSE))  

df <- mmse.merged

# === Add demographics ==========

# age: find the age at earliest baseline scan
# calculate difference from that age

min.ages <- ptdemog %>%
  select(RID, USERDATE, AGE, PTGENDER) %>%
  rename(DateDemogBL=USERDATE, AgeBL=AGE, Gender=PTGENDER) %>%
  mutate(DateDemogBL=as_datetime(ymd(DateDemogBL)),
         AgeBL = as.numeric(AgeBL)) %>%
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

demog.csv <- read.csv('../rawdata/PTDEMOG.csv') %>%
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

demog.csv <- read.csv('../rawdata/PTDEMOG.csv') %>%
  select(RID, PTGENDER) %>%
  filter(RID %in% gender.missing$RID) %>%
  mutate(Gender.Missing=recode(PTGENDER, `1`='Male', `2`='Female')) %>%
  select(-PTGENDER)

df.gender <- left_join(df, demog.csv, by='RID') %>%
  mutate(Gender = ifelse(is.na(Gender), Gender.Missing, Gender)) %>%
  select(-Gender.Missing)

df <- df.gender

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

# === Save ===========

write.csv(df, '../derivatives/main_data_ADNI.csv', quote=F, na='', row.names=F)
