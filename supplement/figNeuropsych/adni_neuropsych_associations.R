# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggsignif))
sh(library(lubridate))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.DATA <- '../../derivatives/adni/data_with_staging.csv'
PATH.ANOVA <- '../../scripts/anova.R'

# === Apply staging to all ADNI data ========

df.all <- read.csv(PATH.DATA)

# select only training data
df <- df.all %>%
  filter(Group == 'TrainingBaseline')

# define colors
stage.colors <- c('0' = 'white', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00',
                  'NS' = 'gray')

# === function for graph/statistics =======

source(PATH.ANOVA)

# === make plots =======

for (y in c('Composite.MEM', 'Composite.EF', 'Composite.LANG')) {
  x = 'PTCStage'
  aov.stats <- my.anova(x, y, df)
  anova.plot(x, y, colors = stage.colors, data = df, sig.y.start = 2)
  ggsave(sprintf('SUPPLEMENT_adni_%s.png', y), width = 6, height = 6, units = 'in')
  write.csv(aov.stats, sprintf('SUPPLEMENT_adni_%s.csv', y))
}

