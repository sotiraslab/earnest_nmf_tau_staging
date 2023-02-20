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

PATH.ADNI <- '../../derivatives/adni/main_data_with_8ptc.csv'
PATH.NMF.8 <- '../../nmf/adni/results/mat/NumBases8.mat'
PATH.NMF.REGIONS <- '../../derivatives/adni/nmf_regions.csv'
PATH.NMF.8.NAMES <- '../../derivatives/adni/names_8ptc.csv'
PATH.ADNI.W.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'

# === Read data ====

df <- read.csv(PATH.ADNI)

# === Get GM Volume in components ==========

mat <- readMat(PATH.NMF.8)
W <- mat$Wnorm

regions <- read.csv(PATH.NMF.REGIONS)$Feature

all.cols <- colnames(df)
vol.cols <- all.cols[grepl('_VOLUME$', all.cols)]

# checks that the order is good
vol.cols.strip = gsub('_VOLUME', '', vol.cols)
regions.strip <- gsub('_SUVR', '', regions)
print(all(regions.strip == vol.cols.strip))

vol.mat <- as.matrix(df[, vol.cols])

project <- vol.mat %*% W

# === Reassign to DF ========

nice.names <-  read.csv(PATH.NMF.8.NAMES)$Component
nice.names <- paste('Cmp.', nice.names, sep='')

colnames(project) <- paste(nice.names, '.Volume', sep='')
df <- cbind(df, project)

# === Get staged order of components =======

stage.order <- read.csv(PATH.ADNI.W.ORDER)$Region

# === Create holder for all heatmap type results =========

stage.order.nice <- gsub('Cmp.', '', stage.order)

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

# === Select X-sectional data ======

df.x <- df %>%
  filter(Group == 'TrainingBaseline')

# === XSect tau vs xsect GM ==========

rs <- init.mat()
ps <- init.mat()

for (i in 1:length(stage.order)) {
  for (j in 1:length(stage.order)) {
    
    suvr <- stage.order[i]
    gm <- paste(stage.order[j], '.Volume', sep='')
    corr <- pcor.test(df.x[[suvr]], df.x[[gm]], list(df.x$Age))
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
        text = element_text(size=15)) +
  scale_fill_colormap(colormap = 'inferno', reverse = T) +
  ylab('Flortaucipir  (SUVR)') +
  xlab(expression(paste('Gray Matter Volume ', (mm^3)))) +
  labs(fill='Pearson R')

ggsave('adni_gm_regression.png', width=8, height=8)