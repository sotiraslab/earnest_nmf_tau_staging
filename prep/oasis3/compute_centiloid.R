# The default OASIS3 table for PUP results has columns for
# centiloid which has a lot of missing values, i.e.
# Centil_fSUVR_rsf_TOT_CORTMEAN.  However, there are columns
# which define the mean cortical PIB/AV45 uptake pre converison,
# e.g. PET_fSUVR_rsf_TOT_CORTMEAN.  This column has no missing values.
# This code just verifies that the non-missing column is giving very
# close values to the provided centiloids.

# There could still be some issue that the missing values are
# due to filtering of QCed data.  Might need to see
# if there is any information on that.

sh <- suppressPackageStartupMessages

sh(library(ggplot2))
sh(library(this.path))

setwd(this.dir())

df = read.csv("../../rawdata/oasis_amyloid.csv")

# https://www.oasis-brains.org/files/OASIS-3_Imaging_Data_Dictionary_v2.0.pdf
centiloid.av45 <- function(x) return((53.6 * x) - 43.2)
centiloid.pib  <- function(x) return((45.0 * x) - 47.5)

df$ManualCentiloid <- ifelse(df$tracer == 'AV45',
                             centiloid.av45(df$PET_fSUVR_rsf_TOT_CORTMEAN),
                             centiloid.pib(df$PET_fSUVR_rsf_TOT_CORTMEAN))

ggplot(data=df, aes(x=Centil_fSUVR_rsf_TOT_CORTMEAN, y=ManualCentiloid)) +
  geom_point(color='blue', alpha=.5) +
  geom_abline(slope=1, intercept = 0)

write.csv(df, '../../derivatives/oasis3/manual_centiloid.csv', quote=F, na='')

