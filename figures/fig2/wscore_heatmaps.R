# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Files needed ========

PATH.ADNI.DATA <- '../../derivatives/adni/main_data_with_8ptc_no_rois.csv'
PATH.PTC.NAMES <- '../../derivatives/adni/names_8ptc.csv'
PATH.OASIS.DATA <- '../../derivatives/oasis3/main_data_with_8ptc_no_rois.csv'
PATH.WSCORE.SCRIPT <- '../../scripts/wscore.R'
PATH.HEATMAP.SCRIPT <- '../../scripts/stage_heatmaps.R'

# === ADNI =========

adni <- read.csv(PATH.ADNI.DATA)

adni.control <- adni[adni$Group == 'ControlBaseline', ]
adni.train <- adni[adni$Group == 'TrainingBaseline', ]

# ADNI - Wscore -------

source(PATH.WSCORE.SCRIPT)

comp.names.short <-  read.csv(PATH.PTC.NAMES)$Component
comp.names <- paste('Cmp.', comp.names.short, sep='')


w.adni <- repeated.wscore.train(control.data = adni.control,
                                y = comp.names,
                                covariates = c('Age', 'Gender'),
                                match.continuous =  c("Age", "CorticalTauAverage"),
                                match.categorical = c("Gender"),
                                portion.train = 0.8,
                                repeats = 200, 
                                seed = T)

# save data with W scores
adni.predicts <- repeated.wscore.predict(w.adni, adni)
colnames(adni.predicts) <- paste(colnames(adni.predicts), '.WScore', sep='')
adni.with.w <- cbind(adni, adni.predicts)

# ADNI - Braak W scores -------
braak.names <- c('BRAAK1_SUVR', 'BRAAK3_SUVR',
                 'BRAAK4_SUVR', 'BRAAK5_SUVR',
                 'BRAAK6_SUVR')

w.adni.braak <- repeated.wscore.train(control.data = adni.control,
                                      y = braak.names,
                                      covariates = c('Age', 'Gender'),
                                      match.continuous =  c("Age", "CorticalTauAverage"),
                                      match.categorical = c("Gender"),
                                      portion.train = 0.8,
                                      repeats = 200, 
                                      seed = T)

adni.braak.predicts <- repeated.wscore.predict(w.adni.braak, adni)
colnames(adni.braak.predicts) <- paste(colnames(adni.braak.predicts), '.W', sep='')

# ADNI - Save -------

save.adni <- cbind(adni, adni.predicts, adni.braak.predicts)
path.out <- '../../derivatives/adni/data_with_wscores.csv'
write.csv(save.adni, path.out, quote=F, na='', row.names=F)

# ADNI - Heatmap -------

source(PATH.HEATMAP.SCRIPT)

thr <- 2.5

wmat <- adni.predicts[adni$Group == 'TrainingBaseline', ]
wmat <- ifelse(wmat > thr, 1, 0)
wdf <- as.data.frame(wmat)
colnames(wdf) <- comp.names.short

df.plot <- cbind(adni.train[, c('RID', 'CDRBinned')], wdf) %>%
  mutate(CDRBinned = recode(CDRBinned,
                            !!! c('0.0' = 'Tau+, CDR=0',
                                  '0.5' = 'Tau+, CDR=0.5',
                                  '1.0+' = 'Tau+, CDR>=1')))

stage.heatmap.by(df.plot, by='CDRBinned',
                 cols=comp.names.short,
                 cats=c('Tau+, CDR=0', 'Tau+, CDR=0.5', 'Tau+, CDR>=1'),
                 colors=c('#440154FF', '#21908CFF', '#FDE725FF'),
                 empty.name = 'Tau-') + 
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust = 1.15,
                                   margin = margin(t = -15, r = 0, b = 0, l = 100)),
        axis.text.y = element_blank(),
        plot.margin = margin(.1,1,.1,1, "cm"),
        plot.title = element_text(margin = margin(t = 0, r = 0, b = -20, l = 0)),
        legend.position = c(.95, .87),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank()) +
  xlab('PTC') +
  ggtitle('ADNI-ADS PTC Positivity')

ggsave('adni_wscore_heatmap.png', width=4, height=10)

# ADNI - Stage Order -------
sumPos <- colSums(wmat)
stage.order <- sort(sumPos, decreasing = T)
names(stage.order) <- gsub('\\.WScore', '', names(stage.order))
stage.df <- data.frame(Region=names(stage.order),
                       NPos=stage.order,
                       row.names = NULL)

path.out <- '../../derivatives/adni/wscore_stage_order.csv'
write.csv(stage.df, path.out)

# === OASIS3 =====

oasis <- read.csv(PATH.OASIS.DATA) %>%
  mutate(Gender=ifelse(GENDER == 1, 'Male', 'Female')) # replaced to match ADNI
oasis.train <- filter(oasis, Group == 'TrainingSet')
oasis.control <- filter(oasis, Group == 'ControlSet')


# OASIS - Wscore -------

source(PATH.WSCORE.SCRIPT)

w.oasis <- repeated.wscore.train(control.data = oasis.control,
                                 y = comp.names,
                                 covariates = c('Age', 'Gender'),
                                 match.continuous =  c("Age", "TotalCtxTauMean"),
                                 match.categorical = c("Gender"),
                                 portion.train = 0.8,
                                 repeats = 200, 
                                 seed = T)

# save data with W scores
oasis.predicts <- repeated.wscore.predict(w.oasis, oasis)
colnames(oasis.predicts) <- paste(colnames(oasis.predicts), '.WScore', sep='')

# OASIS - Braak W scores -------

braak.names <- c('BRAAK1_SUVR', 'BRAAK3_SUVR',
                 'BRAAK4_SUVR', 'BRAAK5_SUVR',
                 'BRAAK6_SUVR')
braak.control.oasis <- oasis.control %>%
  drop_na(all_of(braak.names))
  
w.oasis.braak <- repeated.wscore.train(control.data = braak.control.oasis,
                                       y = braak.names,
                                       covariates = c('Age', 'Gender'),
                                       match.continuous =  c("Age", "TotalCtxTauMean"),
                                       match.categorical = c("Gender"),
                                       portion.train = 0.8,
                                       repeats = 200, 
                                       seed = T)

oasis.braak.predicts <- repeated.wscore.predict(w.oasis.braak, oasis)
colnames(oasis.braak.predicts) <- paste(colnames(oasis.braak.predicts), '.W', sep='')

# OASIS - Save -------
save.oasis <- cbind(oasis, oasis.predicts, oasis.braak.predicts)
path.out <- '../../derivatives/oasis3/data_with_wscores.csv'
write.csv(save.oasis, path.out, quote=F, na='', row.names=F)

# OASIS - Heatmap -------_

source(PATH.HEATMAP.SCRIPT)

wmat <- oasis.predicts[oasis$Group == 'TrainingSet', ]
wmat <- ifelse(wmat > thr, 1, 0)
wdf <- as.data.frame(wmat)
colnames(wdf) <- comp.names.short

df.plot <- cbind(oasis.train[, c('Subject', 'CDR')], wdf) %>%
  mutate(CDR = as.character(CDR),
         CDR = recode(CDR, '0'='Tau+, CDR=0', '0.5'='Tau+, CDR=0.5', '1'='Tau+, CDR>=1'))

stage.heatmap.by(df.plot, by='CDR',
                 cols=comp.names.short,
                 cats=c('Tau+, CDR=0', 'Tau+, CDR=0.5', 'Tau+, CDR>=1'),
                 colors=c('#440154FF', '#21908CFF', '#FDE725FF'),
                 empty.name = 'Tau-') + 
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, vjust = 1.15,
                                   margin = margin(t = -15, r = 0, b = 0, l = 100)),
        axis.text.y = element_blank(),
        plot.margin = margin(.1,1,.1,1, "cm"),
        plot.title = element_text(margin = margin(t = 0, r = 0, b = -20, l = 0)),
        legend.position = c(.95, .87),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank()) +
  xlab('PTC') +
  ggtitle('OASIS3-ADS PTC Positivity')

ggsave('oasis_wscore_heatmap.png', width=4, height=10)

# OASIS - Stage Order -------

sumPos <- colSums(wmat)
stage.order <- sort(sumPos, decreasing = T)
names(stage.order) <- gsub('\\.WScore', '', names(stage.order))
stage.df <- data.frame(Region=names(stage.order),
                       NPos=stage.order,
                       row.names = NULL)

path.out <- '../../derivatives/oasis3/wscore_stage_order.csv'
write.csv(stage.df, path.out)

