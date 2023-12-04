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

# === Apply staging to all ADNI data ========

df.all <- read.csv(PATH.DATA)

# select only training data
df <- df.all %>%
  filter(Group == 'TrainingBaseline')

# define colors
stage.colors <- c('0' = 'white', '1'= '#009E73', '2' = '#F0E442', '3' = '#E69F00', '4' = '#D55E00',
                  'NS' = 'gray')

# === MEM - PTC Staging =========

# run stats, get results and convert to arguments understood by ggsignif
anova <- aov(Composite.LANG ~ PTCStage, data = df)
posthoc <- as.data.frame(TukeyHSD(anova, method='fdr')$PTCStage)

posthoc.res <- posthoc %>%
  rownames_to_column('comparison') %>%
  mutate(annotation = cut(`p adj`,
                          breaks = c(0, 0.001, 0.01, 0.05, Inf),
                          labels = c('***', "**", "*", ""),
                          include.lowest = T)
         )
posthoc.sig <- filter(posthoc.res, `p adj` < 0.05)

comparisons <- str_split(posthoc.sig$comparison, '-')
n.sig <- nrow(posthoc.sig)

mean.pacc <- group_by(df, PTCStage) %>%
  summarise(Composite.LANG = mean(Composite.LANG, na.rm=T)) %>%
  mutate(x=seq(0.5, by=1, length.out=6),
         xend=seq(1.5, by=1, length.out=6))

ggplot(data = df, aes(x = PTCStage, y = Composite.LANG, fill = PTCStage)) +
  geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
  geom_segment(data=mean.pacc, aes(x=x, xend=xend, y=Composite.LANG, yend=Composite.LANG),
               color='black',
               linewidth=1) + 
  scale_fill_manual(values=stage.colors) +
  theme_light() +
  xlab("Stage") +
  theme(legend.position = 'none',
        text = element_text(size=20),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  coord_cartesian(ylim = c(-10, 6)) +
  scale_y_continuous(breaks=c(0, -2.5, -5, -7.5, -10)) +
  geom_signif(comparisons=comparisons,
              annotations = posthoc.sig$annotation,
              y_position = c(1, 2, 3, 4, 5, 1, 2),
              tip_length = 0.01,
              size=.75,
              textsize = 7) +
  ylab('PACC')

# ggsave('adni_pacc_scatter.png', width=6, height=8)

print(summary(anova))

# save posthoc stats
# posthoc.res <- posthoc.res %>%
#   mutate(across(where(is.numeric), round, 3))
# write.csv(posthoc.res, 'SUPPLEMENT_pacc_posthoc_adni.csv')
