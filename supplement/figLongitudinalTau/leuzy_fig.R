# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(colormap))
sh(library(ggsignif))
sh(library(gtools))
sh(library(lme4))
sh(library(lubridate))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required files =========

PATH.DATA <- '../../derivatives/adni/data_with_staging.csv'
PATH.ORDER <- '../../derivatives/adni/wscore_stage_order.csv'

# === read =======

df <- read.csv(PATH.DATA)
ptc.order <- read.csv(PATH.ORDER)$Region

# === compute uptake in each stage ========

df$Stage1 <- df$Cmp.MedialTemporal
df$Stage2 <- rowMeans(df[, c('Cmp.LeftParietalTemporal',
                            'Cmp.RightParietalTemporal',
                            'Cmp.Precuneus')])
df$Stage3 <- rowMeans(df[, c('Cmp.LateralFrontal',
                             'Cmp.Occipital')])
df$Stage4 <- rowMeans(df[, c('Cmp.Sensorimotor',
                             'Cmp.Orbitofrontal')])

stage.cols <- c('Stage1', 'Stage2', 'Stage3', 'Stage4')

# === create disease category assignments =======

df.disease <- df %>%
  group_by(RID) %>%
  mutate(
    DiseaseType = NA,
    DiseaseType = ifelse(first(Group) == 'ControlBaseline', 'A-/CU', DiseaseType),
    DiseaseType = ifelse(first(AmyloidStatus) == 0 & ! first(CDRBinned) == '0.0', 'non-AD', DiseaseType),
    DiseaseType = ifelse(first(AmyloidStatus) == 1 & first(CDRBinned) == '0.0', 'A+/CDR=0', DiseaseType),
    DiseaseType = ifelse(first(AmyloidStatus) == 1 & first(CDRBinned) == '0.5', 'A+/CDR=0.5', DiseaseType),
    DiseaseType = ifelse(first(AmyloidStatus) == 1 & first(CDRBinned) == '1.0+', 'A+/CDR=1.0+', DiseaseType),
  ) %>%
  drop_na(DiseaseType) %>%
  ungroup()

# === select bl data =======

df.bl <- df.disease %>%
  group_by(RID) %>%
  slice(1) %>%
  ungroup()


# ==== longitudinal change in PTC uptake ========

all.cols <- colnames(df.disease)

# for calculating change like Leuzy et al 2022
# for (y in stage.cols) {
#   new.name <- paste(y, '.PercentChange', sep='')
#   long.data <- df.disease %>%
#     select(RID, DateTau, Age, !!sym(y), DiseaseType) %>%
#     group_by(RID) %>%
#     filter(n() >= 2) %>%
#     mutate(DeltaTauDate = as.numeric(difftime(DateTau, first(DateTau), units='days')) / 365.25,
#            BaselineTau = first(!!sym(y)),
#            FollowupTau = !!sym(y),
#            PercentChange = ((FollowupTau - BaselineTau) / BaselineTau) * (100 / DeltaTauDate),
#            PercentChange = ifelse(is.na(PercentChange), 0, PercentChange)) %>%
#     group_by(RID) %>%
#     filter(row_number() == 2)
# 
#   merger <- long.data %>%
#     select(RID, PercentChange)
#   colnames(merger) <- c('RID', new.name)
# 
#   df.bl <- left_join(df.bl, merger, by='RID')
# }

for (y in stage.cols) {
  long.data <- df.disease %>%
    select(RID, DateTau, Age, !!sym(y), DiseaseType) %>%
    group_by(RID) %>%
    filter(n() >= 2) %>%
    mutate(DeltaTauDate = as.numeric(difftime(DateTau, first(DateTau), units='days')) / 365.25,
           LongAge = Age + DeltaTauDate)

  fml <- as.formula(sprintf('%s ~ DeltaTauDate + (1+DeltaTauDate|RID)', y))
  m <- lmer(fml, data = long.data)
  new.name <- paste(y, '.Predicted', sep='')
  long.data[[new.name]] <- predict(m, long.data)

  # plot showing model
  # p <- ggplot(long.data, aes(x=LongAge, y=!!sym(y))) +
  #   geom_point(aes(color=PTCStage), alpha = .7) +
  #   geom_line(aes(y=!!sym(new.name), group=RID, color=PTCStage), alpha= .7)
  # print(p)

  coefs <- coef(m)$RID %>%
    select(DeltaTauDate) %>%
    rownames_to_column(var="RID") %>%
    mutate(RID=as.numeric(RID))
  new.name <- paste(y, '.Slope', sep='')
  colnames(coefs) <- c('RID', new.name)

  df.bl <- left_join(df.bl, coefs, by = 'RID')
}

# variables
stages <- c('A-/CU', 'A+/CDR=0', 'A+/CDR=0.5', 'A+/CDR=1.0+')
plot.columns <- paste(stage.cols, '.Slope', sep='')
colors <- stage.colors <- c('Stage1' = '#009E73',
                            'Stage2' = '#F0E442',
                            'Stage3' = '#E69F00',
                            'Stage4' = '#D55E00')

# # ---- statistics vs stage 0 ------
# 
# compare.stages <- c('2', '3', '4', 'NS')
# compare.ptcs <- plot.columns
# n.tests <- length(compare.stages) * length(compare.ptcs)
# t.data.vs.0 <- data.frame(StageA = rep('0', n.tests),
#                      StageB = rep(NA, n.tests),
#                      PTC = rep(NA, n.tests),
#                      tval = rep(NA, n.tests),
#                      pval = rep(NA, n.tests))
# 
# for (i in seq_along(compare.stages)) {
#   stage <- compare.stages[i]
#   x.data <- filter(df.bl, PTCStage == '0')
#   y.data <- filter(df.bl, PTCStage == stage)
#   
#   for (j in seq_along(compare.ptcs)) {
#     ptc <- compare.ptcs[j]
#     x <- x.data[[ptc]]
#     y <- y.data[[ptc]]
#     comparison <- t.test(x, y, alternative = 'less')
#     
#     row <- (i-1) * length(compare.ptcs) + j
#     t.data.vs.0[row, 'StageB'] <- stage
#     t.data.vs.0[row, 'PTC'] <- ptc
#     t.data.vs.0[row, 'tval'] <- comparison$statistic
#     t.data.vs.0[row, 'pval'] <- comparison$p.value
#   }
# }
# 
# # ---- statistics vs stage previous ------
# 
# compare.stages <- c('1', '2', '3', '4')
# compare.ptcs <- plot.columns
# n.tests <- length(compare.stages) * length(compare.ptcs)
# t.data.vs.prev <- data.frame(StageA = rep(NA, n.tests),
#                      StageB = rep(NA, n.tests),
#                      PTC = rep(NA, n.tests),
#                      tval = rep(NA, n.tests),
#                      pval = rep(NA, n.tests))
# 
# for (i in seq_along(compare.stages)) {
#   stageB <- compare.stages[i]
#   stageA <- as.character(as.integer(stageB) - 1)
#   x.data <- filter(df.bl, PTCStage == stageA)
#   y.data <- filter(df.bl, PTCStage == stageB)
#   
#   for (j in seq_along(compare.ptcs)) {
#     ptc <- compare.ptcs[j]
#     x <- x.data[[ptc]]
#     y <- y.data[[ptc]]
#     comparison <- t.test(x, y, alternative = 'less')
#     
#     row <- (i-1) * length(compare.ptcs) + j
#     t.data.vs.prev[row, 'StageA'] <- stageA
#     t.data.vs.prev[row, 'StageB'] <- stageB
#     t.data.vs.prev[row, 'PTC'] <- ptc
#     t.data.vs.prev[row, 'tval'] <- comparison$statistic
#     t.data.vs.prev[row, 'pval'] <- comparison$p.value
#   }
# }
# 
# 
# # ---- combine T data -------
# 
# t.data.all <- rbind(t.data.vs.prev, t.data.vs.0) %>%
#     mutate(pval.correct = p.adjust(pval, method = 'fdr'),
#            annotation = cut(
#              pval.correct,
#              breaks = c(0, 0.001, 0.01, 0.05, Inf),
#              labels = c('***', "**", "*", ""),
#              include.lowest = T),
#            PTC = str_replace_all(PTC, 'Cmp.|.Slope', ''),
#            ) 
# 
# write.csv(t.data.all, 'longitudinal_tau_t_stats.csv',
#           quote = F, na = '', row.names = F)

# ==== plot ========

# panel for each stage
# bar for each PTC

# plot loop
for (stage in stages) {
  n.stage <- sum(df.bl$DiseaseType == stage)
  plot.data <- df.bl %>%
    filter(DiseaseType == stage) %>%
    select(all_of(plot.columns)) %>%
    pivot_longer(everything(), names_to = 'PTCStage', values_to = 'TauChange') %>%
    mutate(PTCStage = str_replace_all(PTCStage, '.Slope', ''),
           PTCStage = factor(PTCStage, levels = stage.cols))

  # get means by group
  n.categories <- length(unique(plot.data$PTCStage))
  means <- group_by(plot.data, PTCStage) %>%
    summarise(Mean = mean(TauChange, na.rm=T)) %>%
    mutate(x=seq(0.5, by=1, length.out=n.categories),
           xend=seq(1.5, by=1, length.out=n.categories))
  
  
  # draw
  p <- ggplot(data = plot.data, aes(x = PTCStage, y = TauChange, fill = PTCStage)) +
    scale_fill_manual(values = colors) +
    geom_point(position = position_jitter(width = 0.2, seed=42), shape=21, size=3) +
    geom_segment(data=means, aes(x=x, xend=xend, y=Mean, yend=Mean),
                 color='black', linewidth=1) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 14),
          legend.position = 'none') +
    ylab('Tau Change (SUVR/y)') +
    ggtitle(sprintf('PTC Stage %s (n=%s)', stage, n.stage)) +
    coord_cartesian(ylim = c(-0.1, .5), expand = F)
  
  # add stats
  # if (stage != '0') {
  #   stats <- t.data.all %>%
  #     filter(StageB == stage & StageA == '0') 
  #   means$annotation <- stats$annotation
  #   means$y <- -0.06
  #   p <- p +
  #     geom_text(data = means,
  #               aes(x = PTC, y = y, label = annotation),
  #               size = 10,
  #               color = 'blue')
  # }
  # 
  # if (stage  %in% c('2', '3', '4')) {
  #   stats <- t.data.all %>%
  #     filter(StageB == stage & StageA == as.character(as.integer(stage) - 1)) 
  #   means$annotation <- stats$annotation
  #   means$y <- -0.08
  #   p <- p +
  #     geom_text(data = means,
  #               aes(x = PTC, y = y, label = annotation),
  #               size = 10,
  #               color = 'red')
  # }

  print(p)
  
  # ggsave(sprintf('tau_change_stage_%s.png', stage), width = 6, height = 6)
  
}





