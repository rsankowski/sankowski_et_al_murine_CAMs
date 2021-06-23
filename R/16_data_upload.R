library(tools)
library(tidyverse)

files1 <- list.files("data", pattern = ".RData")
map(files1, function(i) {
  a <- md5sum(file.path("data", i))
}) %>% 
  unlist %>%
  data.frame() %>%
  write.csv(file.path("data", "checksums_rdata.csv"))

files2 <- list.files(file.path("data", "counts"), pattern = ".coutt")
map(files2, function(i) {
  a <- md5sum(file.path("data", "counts", i))
}) %>% 
  unlist %>%
  data.frame() %>%
  write.csv(file.path("data", "checksums_counts.csv"))
