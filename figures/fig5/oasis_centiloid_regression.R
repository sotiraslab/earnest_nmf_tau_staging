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

ggplot(data=long.centiloid, aes(x=Centiloid, y=AvgSUVR, group=PTC, color=PTC, linetype=PTC)) + 
  geom_smooth(alpha=.3, method='lm') +
  scale_color_manual(values=colors) +
  ylab('Flortaucipir (Z-score)') +
  theme_light() +
  theme(text = element_text(size=20),
        legend.key.size = unit(2, 'line')) +
  coord_cartesian(xlim=c(15, 200), ylim=c(0, 10)) +
  scale_x_continuous(breaks=c(25, 50, 75, 100, 125, 150, 175, 200))

ggsave(filename = 'centiloid_loadings_oasis.png', width=12, height=6)

# === stats ========

sh(library(emmeans))

m <- lm(AvgSUVR ~ Centiloid*PTC + Age + GENDER, data=long.centiloid)
summary(m)
em <- emmeans(m, 'PTC')
em.summary <- summary(pairs(em, adjust='fdr')) %>%
  mutate(across(where(is.numeric), round, 3),
         annotation = cut(p.value,
                          breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c('***', "**", "*", ""),
                          include.lowest = T)) %>%
  select(-df, -SE)

# save
write.csv(em.summary, 'SUPPLEMENT_emmeans_centiloid_oasis3.csv', row.names = F)

# === Supplement: Stats figure ==========

sh(library(stringr))

split.comps <- as.data.frame(str_split(em.summary$contrast, ' - ', simplify = T))
colnames(split.comps) <- c('RegionA', 'RegionB')

plot.data <- cbind(split.comps, em.summary) %>%
  mutate(RegionA = factor(RegionA, levels=rev(ordered.nice.names)),
         RegionB = factor(RegionB, levels=ordered.nice.names))

ggplot() +
  geom_tile(data=plot.data, mapping=aes(y=RegionA, x=RegionB, fill=t.ratio)) +
  coord_equal() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill='white'),
        text = element_text(size=20),
        axis.ticks = element_blank()) +
  geom_text(data=plot.data, mapping=aes(y=RegionA, x=RegionB, label=annotation),
            color='white', size=7) +
  xlab('PTC') +
  ylab('PTC')

ggsave('SUPPLEMENT_centiloid_regression_stats_oasis.png', width=8, height=8)