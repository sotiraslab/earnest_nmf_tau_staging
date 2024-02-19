
# This script tests that the functions in ptc_staging.R are
# able to reproduce the results presented in the manuscript.

# This can be run so long as `run_build_datasets.sh` has been run.

library(dplyr)
library(this.path)


# Read data -----
setwd(this.dir())

path.adni <- '../derivatives/adni/main_data.csv'
path.oasis <- '../derivatives/oasis3/main_data.csv'

adni <- read.csv(path.adni) %>% filter(Group == 'TrainingBaseline')
oasis <- read.csv(path.oasis) %>% filter(Group == 'TrainingSet')

tau.rois <- colnames(adni)[grepl('CTX_.*_SUVR', colnames(adni), perl = T)]

# Load functions -----

source('ptc_staging.R')

# 1. ADNI with pipeline -----
adni$Sex <- ifelse(adni$Gender == 'Male', 1, 0)
adni.results <- ptc.staging.pipeline(adni,
                                     tau.roi.columns = tau.rois)

# 2. ADNI "by hand" -----
adni.uptakes <- ptc.uptake(adni, tau.rois)
adni.uptakes$Age <- adni$Age
adni.uptakes$Sex <- adni$Sex
adni.positivity <- ptc.wscores(adni.uptakes, cutoff = 2.5)
adni.stages <- ptc.staging(adni.positivity)

print('ADNI PTC positivity:')
print(colSums(adni.positivity))

print('ADNI stages:')
print(table(adni.stages))

# 3. OASIS with pipeline -----
oasis$Sex <- ifelse(oasis$GENDER == 1, 1, 0)
oasis.results <- ptc.staging.pipeline(oasis,
                                     tau.roi.columns = tau.rois,
                                     model = 'oasis')

# 4. OASIS "by hand" -----
oasis.uptakes <- ptc.uptake(oasis, tau.rois)
oasis.uptakes$Age <- oasis$Age
oasis.uptakes$Sex <- oasis$Sex
oasis.positivity <- ptc.wscores(oasis.uptakes, cutoff = 2.5, model = 'oasis')
oasis.stages <- ptc.staging(oasis.positivity)

print('OASIS PTC positivity:')
print(colSums(oasis.positivity))

print('OASIS stages:')
print(table(oasis.stages))

# Scratch for the readme -----

source('ptc_staging.R')
path <- '../derivatives/adni/main_data.csv'
df <- read.csv(path)
df <- df[df$Group == 'TrainingBaseline', ]

all.columns <- colnames(df)
rois <- all.columns[grepl('CTX_.*_SUVR', all.columns)]
tau.uptakes <- ptc.uptake(df, tau.roi.columns = rois)

# adding age & sex to the uptakes table
tau.uptakes$Age <- df$Age
tau.uptakes$Sex <- ifelse(df$Gender == 'Male', 1, 0)

tau.wscores <- ptc.wscores(tau.uptakes,
                           model = 'adni',
                           cutoff = 2.5,
                           age.column = 'Age',
                           sex.column = 'Sex')

tau.stages <- ptc.staging(tau.wscores)

tau.uptakes <- ptc.uptake(df, tau.roi.columns = rois)
my.binary <- as.data.frame(ifelse(tau.uptakes >= 1.20, 1, 0))
my.stages <- ptc.staging(my.binary)
table(my.stages)
