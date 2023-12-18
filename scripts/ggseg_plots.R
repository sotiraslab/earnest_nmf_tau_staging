
library(colormap)
library(dplyr)
library(ggseg)
library(gsubfn)
library(tidyverse)

check.in.ggseg <- function(label) {
  return(label %in% brain_labels(dk) | label %in% brain_labels(aseg))
}

adni.labels.to.ggseg <- function(labels, filter.failures = T,
                                 warn = T) {
  
  # cortical
  cortical <- labels %>%
    str_replace('^CTX_', '') %>%
    str_replace('_VOLUME$', '') %>%
    str_replace('_SUVR$', '') %>%
    tolower()
  
  # subcortical
  subcortical <- labels %>%
    str_replace('_SUVR$', '') %>%
    str_replace('_VOLUME$', '') %>%
    str_replace_all('_', '-') %>%
    tolower() %>%
    str_replace('dc$', 'DC')
  subcortical <- gsubfn('^.|-.', toupper, subcortical)
  
  # assemble output
  output <- ifelse(str_detect(labels, '^CTX'), cortical, subcortical)
  
  # check success
  failures <- ! check.in.ggseg(output)
  
  # filter
  if (filter.failures) {
    output <- ifelse(failures, labels, output)
  }
  
  # warning
  if (warn) {
    issues <- output[failures]
    if (length(issues) > 0) {
      s <- paste(issues, collapse=', ')
      print(sprintf('WARNING! labels not matching ggseg: %s', s))
    }
  }
  
  return(output)
}

plot.cortex<- function(values, regions, vmin = NULL, vmax = NULL,
                       name = 'My Plot', legend = T, cm = 'inferno') {
  if (is.null(vmin)) vmin <- min(values)
  if (is.null(vmax)) vmax <- max(values)
  
  df <- data.frame(label=regions, value=values)
  
  # cortical
  atlas <- data.frame(label = brain_labels(dk)) %>%
    left_join(df, by='label') %>%
    drop_na()
  
  plot.cortical <- ggplot(atlas) +
    geom_brain(atlas = dk,
               color='black',
               size=.5,
               aes(fill = value)) +
    theme_void() +
    ggtitle(name) + 
    scale_fill_colormap(colormap=cm, limits=c(vmin, vmax), oob = scales::squish)
  
  if (! legend) plot.cortical <- plot.cortical + theme(legend.position = 'none')
  
  return(plot.cortical)
}