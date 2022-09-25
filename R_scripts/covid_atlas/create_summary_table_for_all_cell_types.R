library(tidyverse)
library(glue)

in_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/summary_table/"
out_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/final_sum_table/"


files = list.files(path = glue("{in_dir}"))
sum_table <- files %>% 
    map(~ readRDS(file.path(in_dir, .))) %>%
    purrr::reduce(rbind)

sum_table %>% saveRDS(glue("{out_dir}merged_sum_table.rds"))

