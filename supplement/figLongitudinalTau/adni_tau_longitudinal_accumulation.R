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

df <- read.csv(PATH.DATA) %>%
  filter(grepl('Training', Group))
df.bl <- df %>%
  filter(Group == 'TrainingBaseline')

ptc.order <- read.csv(PATH.ORDER)$Region

# ==== longitudinal change in PTC uptake ========

all.cols <- colnames(df)

# for (y in ptc.order) {
#   new.name <- paste(y, '.PercentChange', sep='')
#   long.data <- df %>%
#     select(RID, DateTau, Age, !!sym(y), PTCStage) %>%
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

for (y in ptc.order) {
  long.data <- df %>%
    select(RID, DateTau, Age, !!sym(y), PTCStage) %>%
    group_by(RID) %>%
    filter(n() >= 2) %>%
    mutate(DeltaTauDate = as.numeric(difftime(DateTau, first(DateTau), units='days')) / 365.25,
           LongAge = Age + DeltaTauDate)

  fml <- as.formula(sprintf('%s ~ DeltaTauDate + (1+DeltaTauDate|RID)', y))
  m <- lmer(fml, data = long.data)
  new.name <- paste(y, '.Predicted', sep='')
  long.data[[new.name]] <- predict(m, long.data)

  p <- ggplot(long.data, aes(x=LongAge, y=!!sym(y))) +
    geom_point(aes(color=PTCStage), alpha = .7) +
    geom_line(aes(y=!!sym(new.name), group=RID, color=PTCStage), alpha= .7)
  print(p)

  coefs <- coef(m)$RID %>%
    select(DeltaTauDate) %>%
    rownames_to_column(var="RID") %>%
    mutate(RID=as.numeric(RID))
  new.name <- paste(y, '.Slope', sep='')
  colnames(coefs) <- c('RID', new.name)

  df.bl <- left_join(df.bl, coefs, by = 'RID')
}

# ==== plot ========

df.bl <- filter(df.bl, ! is.na(Cmp.MedialTemporal.Slope))
plot.cols <- paste(ptc.order, '.Slope', sep='')


for (y in plot.cols) {
  p <- ggplot(df.bl, aes(x=PTCStage, y=!!sym(y))) + 
    geom_bar(stat = 'summary', fun = mean)
  print(p)
}



