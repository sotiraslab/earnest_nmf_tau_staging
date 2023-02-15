# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggplot2))
sh(library(tidyverse))
sh(library(this.path))

# === Required files =========

PATH.OASIS.ZSCORE <- '../../derivatives/oasis3/data_loadings_zscore.csv'
PATH.ADNI.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'

# === Load Data ========

df <- read.csv(PATH.OASIS.ZSCORE) %>%
  filter(Group == 'TrainingSet')

cols <- colnames(df)
components <- cols[grepl('Cmp', cols)]

nice.names <- gsub('Cmp.', '', components)
renamer <- nice.names
names(renamer) <- components

# read component order  ========

order.comps <- read.csv(PATH.ADNI.ORDER)$Region
ordered.nice.names <- gsub('Cmp.', '', order.comps)

# === plot ========

long.centiloid <- df %>%
  select(all_of(components), "Centiloid", Age, GENDER) %>%
  pivot_longer(all_of(components), values_to='AvgSUVR', names_to='PTC') %>%
  mutate(PTC=recode(PTC, !!!renamer),
         PTC=factor(PTC, levels=ordered.nice.names))

colors = colormap('jet', nshades = 8)
colors[1] <- 'black'
names(colors) <- ordered.nice.names

ggplot(data=long.centiloid, aes(x=Centiloid, y=AvgSUVR, group=PTC, color=PTC)) + 
  geom_smooth(alpha=.3, method='lm') +
  scale_color_manual(values=colors) +
  ylab('Flortaucipir (Z-score)') +
  theme_light() +
  theme(text = element_text(size=20)) +
  coord_cartesian(xlim=c(15, 200), ylim=c(0, 10)) +
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150, 175, 200))

ggsave(filename = 'centiloid_loadings_oasis.png', width=12, height=6)

# === stats ========

sh(library(emmeans))

m <- lm(AvgSUVR ~ Centiloid + PTC + Age + GENDER, data=long.centiloid)
em <- emmeans(m, 'PTC')
em.summary <- summary(pairs(em, adjust='fdr'))

# save
write.csv(em.summary, 'emmeans_centiloid_oasis.csv', row.names = F)
