# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggplot2))
sh(library(ggseg))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# ==== Required files ========

PATH.ADNI.ORDER <- '../fig3/adni_wscore_stage_order.csv'
PATH.REGION.ASSIGNMENT <- '../fig3/assignments.csv'

# ==== Plot adni ====

adni.pos <- read.csv(PATH.ADNI.ORDER)
assignments <- read.csv(PATH.REGION.ASSIGNMENT)

# manually define the staging system
stages <- c('Cmp.MedialTemporal'=1,
            'Cmp.LeftParietalTemporal'=2,
            'Cmp.RightParietalTemporal'=2,
            'Cmp.Precuneus'=2,
            'Cmp.Occipital'=3,
            'Cmp.LateralFrontal'=3,
            'Cmp.Orbitofrontal'=4,
            'Cmp.Sensorimotor'=4)
merger <- data.frame(name=names(stages),
                     Stage=factor(unname(stages)))

plot.data <- left_join(assignments, merger, by='name')
colors = c('1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00')

p <- ggplot(plot.data) +
  geom_brain(atlas = dk,
             aes(fill = Stage),
             position = position_brain(hemi ~ side ),
             color='black',
             linewidth=4) +
  theme_void() +
  theme(text = element_text(size=20)) +
  scale_fill_manual(values=colors)
p

ggsave('stage_map.png', p, width=7, height=7)


