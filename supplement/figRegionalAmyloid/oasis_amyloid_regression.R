# === Imports ============

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ppcor))
sh(library(R.matlab))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set Working Directory ============

setwd(this.dir())

# === Required Files ============

PATH.DATA <- '../../derivatives/oasis3/main_data_with_8ptc.csv'
PATH.OASIS.AMYLOID <- '../../rawdata/oasis_amyloid.csv'
PATH.NMF.8 <- '../../nmf/adni/results/mat/NumBases8.mat'
PATH.NMF.REGIONS <- '../../derivatives/adni/nmf_regions.csv'
PATH.PTC.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'

# === Read data ====

df <- read.csv(PATH.DATA) %>%
  filter(Group == 'TrainingSet')

regions <- read.csv(PATH.NMF.REGIONS)$Feature
ptcs <- read.csv(PATH.PTC.ORDER)$Region

# === preprocess ADNI amyloid data ========

# extract only the subjects
adrc_session_to_number <- function(col) {
  extracted <- str_extract(col, 'd\\d+')
  no.leading.d <- substr(extracted, 2, nchar(extracted))
  number <- as.numeric(no.leading.d)
  
  return (number)
}

oasis.amy <- read.csv(PATH.OASIS.AMYLOID)
all.cols <- colnames(oasis.amy)
roi.cols <- all.cols[str_detect(all.cols, 'PET_fSUVR_(L|R)_.*') &
                       ! str_detect(all.cols, 'CRPCLM')]
oasis.amy <- oasis.amy %>%
  mutate(
    Subject = str_extract(PUP_PUPTIMECOURSEDATA.ID, 'OAS\\d+'),
    SessionAmyloid = adrc_session_to_number(PUP_PUPTIMECOURSEDATA.ID),
    AmyloidTracer = tracer,
         ) %>%
  select(Subject, SessionAmyloid, AmyloidTracer, all_of(roi.cols))

# verify roi order is the same
converter <- data.frame(ADNI = regions,
                        OASIS = roi.cols)
colnames(oasis.amy) <- c('Subject', 'SessionAmyloid', 'AmyloidTracer', regions)

# === Merge into dataset ========

abeta.big <- df %>%
  select(Subject, SessionAmyloid, AmyloidTracer, Age, all_of(ptcs)) %>%
  left_join(oasis.amy, by = c('Subject', 'SessionAmyloid', 'AmyloidTracer'))

# === Get amyloid SUVR in components ==========

mat <- readMat(PATH.NMF.8)
W <- mat$Wnorm
amy.mat <- as.matrix(abeta.big[, regions])
project <- amy.mat %*% W

# === Reassign to DF ========

amy.ptcs <- paste(ptcs, '.Amyloid', sep='')
colnames(project) <- amy.ptcs
abeta <- abeta.big %>%
  select(-all_of(regions)) %>%
  cbind(project)

# === Create holder for all heatmap type results =========

stage.order.nice <- gsub('Cmp.', '', ptcs)

init.mat <- function() {
  mat <- as.data.frame(matrix(nrow=8, ncol=8))
  rownames(mat) <- stage.order.nice
  colnames(mat) <- stage.order.nice
  return(mat)
}

melt.mat <- function(mat) {
  mat <- mat %>%
    rownames_to_column('component1') %>%
    pivot_longer(all_of(colnames(mat)), names_to = 'component2', values_to = "value") %>%
    mutate(component1=factor(component1, levels=rev(stage.order.nice)),
           component2=factor(component2, levels=stage.order.nice))
  
  return(mat)
}

# === wrap method into function ========


amy.tau.heatmap <- function(data, saveplot=NULL, savestats=NULL,
                            xlab = 'Amyloid-PET (SUVR)') {
  rs <- init.mat()
  ps <- init.mat()
  
  for (i in 1:length(ptcs)) {
    for (j in 1:length(ptcs)) {
      
      suvr <- ptcs[i]
      amy <- paste(ptcs[j], '.Amyloid', sep='')
      corr <- pcor.test(data[[suvr]], data[[amy]], list(data$Age))
      rs[i, j] <- corr$estimate
      ps[i, j] <- corr$p.value
      
    }
  }
  
  bup <- rs
  rs <- melt.mat(rs)
  ps <- melt.mat(ps)
  ps$value <- p.adjust(ps$value, method='fdr')
  
  rs.plot <- rs[ps$value <= 0.05, ]
  
  ggplot() + 
    geom_tile(data=rs.plot, aes(x=component2, y=component1, fill=value)) + 
    coord_equal() +
    theme_classic() + 
    theme(axis.text.x = element_text(angle=45, hjust=1),
          text = element_text(size=15),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_colormap(colormap = 'bathymetry') +
    ylab('Flortaucipir  (SUVR)') +
    xlab(xlab) +
    labs(fill='Pearson R')
  
  if (! is.null(saveplot)) {
    ggsave(saveplot, width=8, height=8)
  }
  
  if (! is.null(savestats)) {
    stats.df <- rs
    colnames(stats.df) <- c('regionFTP', 'regionAmy', 'R')
    stats.df <- stats.df %>%
      mutate(p.value = round(ps$value, 7),
             annotation = cut(p.value,
                              breaks = c(0, 0.001, 0.01, 0.05, Inf),
                              labels = c('***', "**", "*", ""),
                              include.lowest = T),
             R = round(R, 3))
    
    write.csv(stats.df, savestats)
    
  }
}

# ==============

amy.tau.heatmap(abeta,
                saveplot = 'oasis_all_data_heatmap.png',
                savestats = 'oasis_all_data_stats.csv',
                xlab = 'Amyloid-PET (SUVR)')

amy.tau.heatmap(filter(abeta, AmyloidTracer == 'AV45'),
                saveplot = 'oasis_av45_heatmap.png',
                savestats = 'oasis_av45_stats.csv',
                xlab = 'Florbetapir (SUVR)')

amy.tau.heatmap(filter(abeta, AmyloidTracer == 'PIB'),
                saveplot = 'oasis_pib_heatmap.png',
                savestats = 'oasis_pib_stats.csv',
                xlab = 'PIB (SUVR)')

# ==========

# x <- 'Cmp.MedialTemporal'
# 
# plot(abeta$Cmp.RightParietalTemporal.Amyloid, abeta[[x]])
# abline(lm(abeta[[x]] ~ abeta$Cmp.RightParietalTemporal.Amyloid))
# 
# plot(abeta$Cmp.LeftParietalTemporal.Amyloid, abeta[[x]])
# abline(lm(abeta[[x]] ~ abeta$Cmp.LeftParietalTemporal.Amyloid))
