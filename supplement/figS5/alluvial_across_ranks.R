# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(anticlust))
sh(library(colormap))
sh(library(ggalluvial))
sh(library(ggplot2))
sh(library(ggseg))
sh(library(gtools))
sh(library(R.matlab))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

# === Set WD ===========

setwd(this.dir())

# === Required Files =========

PATH.MAT <- '../../nmf/adni/results/mat/'
PATH.REGIONS <- '../../derivatives/adni/nmf_regions_ggseg.csv'

# ==== Read data =========

# NMF solutions
mats <- mixedsort(list.files(PATH.MAT, full.names = T))
ranks <- as.numeric(str_extract(mats, '\\d+'))
  
# NMF regions
regions <- read.csv(PATH.REGIONS)

# holder for NMF winner take all assignments
assignments.mat <- matrix(data = NA, nrow = nrow(regions), ncol = length(ranks))

# loop to read and do winner take all
for (i in seq_along(mats)) {
  path <- mats[i]
  print(path)
  mat.file <- readMat(path)
  Wnorm <- mat.file$Wnorm
  wta <- apply(Wnorm, 1, which.max)
  assignments.mat[, i] <- wta
}

assignments <- as.data.frame(assignments.mat)
rank.cols <- as.character(ranks)
colnames(assignments) <- rank.cols
assignments$Region <- regions$label


# ==== Create plot data ===========

plot.alluvial <- function(ranks) {
  plot.data <- assignments %>%
    pivot_longer(all_of(rank.cols), names_to = 'Rank', values_to = 'Component') %>%
    mutate(Rank = as.numeric(Rank),
           Component = factor(Component)) %>%
    filter(Rank %in% ranks)
  
  
  ggplot(plot.data, aes(x = Rank, stratum = Component, alluvium = Region,
                        fill = Component, label = Component)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", aes.flow='forward') +
    geom_stratum() +
    theme_light()
}

ranks <- c(2, 6, 8)
plot.data <- assignments %>%
  pivot_longer(all_of(rank.cols), names_to = 'Rank', values_to = 'Component') %>%
  mutate(Rank = as.numeric(Rank),
         Component = factor(Component)) %>%
  filter(Rank %in% ranks)

t <- plot.data %>%
  filter(Group == 8) %>%
  arrange(Component)

plot.data$Region <- factor(plot.data$Region, levels=t$Region)

ggplot(plot.data, aes(x = Rank, stratum = Component, alluvium = Region,
                      fill = Component, label = Component)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes.flow='forward') +
  geom_stratum() +
  theme_light()

# plot.alluvial(c(2, 4, 6, 8, 10))

