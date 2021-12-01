# mengyao
library(tidyverse)
library(tidyseurat)
devtools::install_github("stemangiola/tidySingleCellExperiment")

main_S %>%
  distinct(sample)


'''
dir_names = c("GSM4138110/data/")
file_names_raw = sprintf("%s", dir_names)

counts = 
  file_names_raw %>%
  DropletUtils::read10xCounts(col.names=T)
counts
counts %>% saveRDS("project/GSM4138110/RDS")
'''
read.scv