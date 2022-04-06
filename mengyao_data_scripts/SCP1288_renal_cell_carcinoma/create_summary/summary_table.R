library(tidyverse)
library(glue)
in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/summary_table/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/final_sum_table/"


files = list.files(path = glue("{in_dir}"))
sum_table <- files %>% 
    map(~ readRDS(file.path(in_dir, .))) %>%
    reduce(rbind)

sum_table %>% saveRDS(glue("{out_dir}sum_table.rds"))
