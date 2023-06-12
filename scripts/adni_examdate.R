library(dplyr)

get.examdate.from.registry <- function(data) {
  
  # select columns which uniquely identify visits
  a <- select(data, RID, COLPROT, VISCODE)
  
  # select same columns + EXAMDATE in registry
  # groupby accounts for some duplicate visits
  b <- select(registry, RID, COLPROT, VISCODE, EXAMDATE) %>%
    group_by(RID, COLPROT, VISCODE) %>%
    slice(n()) %>%
    ungroup()
  
  # left join to get date
  c <- left_join(a, b, by=c('RID', 'COLPROT', 'VISCODE'))
  
  return (as.character(c$EXAMDATE))
}