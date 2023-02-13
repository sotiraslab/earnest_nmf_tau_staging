# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Necessary files ========

PATH.ADNI <- '../../derivatives/adni/main_data.csv'
PATH.OASIS3 <- '../../derivatives/oasis3/main_data.csv'

adni <- read.csv(PATH.ADNI)
adni.ads <- filter(adni, Group == 'TrainingBaseline')
adni.cn <- filter(adni, Group == 'ControlBaseline')

oasis <- read.csv(PATH.OASIS3)
oasis.ads <- filter(oasis, Group == 'TrainingSet')
oasis.cn <- filter(oasis, Group == 'ControlSet')

# === Descriptives =========

rows <- c('N', 'Age.mean', 'Age.sd', 'Sex.M', 'Sex.F',
           'CDR.0', 'CDR.0.5', 'CDR.1.0+', 'MMSE.mean',
           'MMSE.sd', 'Centiloid.mean', 'Centiloid.sd',
           'APOE.E4+', 'APOE.E4-', 'APOE.NA')
nas <- rep(NA, length(rows))

df <- data.frame(ADNI.ADS = nas, OASIS3.ADS = nas, row.names = rows)

# === fill in ADNI ==========

col <- 'ADNI.ADS'
x <- adni.ads

df['N', col] <- nrow(x)
df['Age.mean', col] <- mean(x$Age)
df['Age.sd', col] <- sd(x$Age)
df['Sex.M', col] <- sum(x$Gender == 'Male')
df['Sex.F', col] <- sum(x$Gender == 'Female')
df['CDR.0', col] <- sum(x$CDRBinned == '0.0')
df['CDR.0.5', col] <- sum(x$CDRBinned == '0.5')
df['CDR.1.0+', col] <- sum(x$CDRBinned == '1.0+')
df['MMSE.mean', col] <- mean(x$MMSE)
df['MMSE.sd', col] <- sd(x$MMSE)
df['Centiloid.mean', col] <- mean(x$Centiloid)
df['Centiloid.sd', col] <- sd(x$Centiloid)
df['APOE.E4+', col] <- sum(x$HasE4, na.rm = T)
df['APOE.E4-', col] <- sum(! x$HasE4, na.rm = T)
df['APOE.NA', col] <- sum(is.na(x$HasE4), na.rm = T)


# === fill in OASIS ===========

col <- 'OASIS3.ADS'
x <- oasis.ads

df['N', col] <- nrow(x)
df['Age.mean', col] <- mean(x$Age)
df['Age.sd', col] <- sd(x$Age)
df['Sex.M', col] <- sum(x$GENDER == 1)
df['Sex.F', col] <- sum(x$GENDER == 2)
df['CDR.0', col] <- sum(x$CDR == 0)
df['CDR.0.5', col] <- sum(x$CDR == .5)
df['CDR.1.0+', col] <- sum(x$CDR == 1)
df['MMSE.mean', col] <- mean(x$MMSE)
df['MMSE.sd', col] <- sd(x$MMSE)
df['Centiloid.mean', col] <- mean(x$Centiloid)
df['Centiloid.sd', col] <- sd(x$Centiloid)
df['APOE.E4+', col] <- sum(x$HasE4, na.rm = T)
df['APOE.E4-', col] <- sum(! x$HasE4, na.rm = T)
df['APOE.NA', col] <- sum(is.na(x$HasE4), na.rm = T)

# === T-Tests ============

x <- adni.ads
y <- oasis.ads

ttests <- list(
  'Age' = t.test(x$Age, y$Age, alternative = 't'),
  'MMSE' = t.test(x$MMSE, y$MMSE, alternative = 't'),
  'Centiloid' = t.test(x$Centiloid, y$Centiloid, alternative = 't')
)

ttest.df <- sapply(ttests, function(x) {
  c('t'=unname(x$statistic),
    'df'=unname(x$parameter),
    'p'=round(unname(x$p.value), 3))
})

ttest.df <- as.data.frame(t(ttest.df))

# === Chi squared =========

x <- adni.ads
y <- oasis.ads

my.chi <- function(v1, v2, label1, label2) {
  df <- data.frame(x=c(v1, v2),
                   label=c(rep(label1, length(v1)), rep(label2, length(v2))))
  return(chisq.test(table(df$label, df$x)))
}

y$Gender <- ifelse(y$GENDER == 1, 'Male', 'Female')
y$CDRBinned <- recode(y$CDR, `0`='0.0', `0.5`='0.5', `1`='1.0+')

chis <- list(
  'Sex' = my.chi(x$Gender, y$Gender, 'ADNI.ADS', 'OASIS3.ADS'),
  'CDR' = my.chi(x$CDRBinned, y$CDRBinned, 'ADNI.ADS', 'OASIS3.ADS'),
  'E4' = my.chi(x$HasE4, y$HasE4, 'ADNI.ADS', 'OASIS3.ADS')
)

chi.df <- sapply(chis, function(x) {
  c('chi'=unname(x$statistic),
    'df'=unname(x$parameter),
    'p'=round(unname(x$p.value), 3))
})

chi.df <- as.data.frame(t(chi.df))

# === Save ========

write.csv(df, 'descriptives.csv')
write.csv(chi.df, 'chi.csv')
write.csv(ttest.df, 'ttest.csv')


