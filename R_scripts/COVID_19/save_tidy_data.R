# modify into tidy_seurat object
library(tidyverse)
library(glue)

in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/all_raw_data/"

in_data = readRDS(glue("{in_dir}s41587-020-0602-4_COVID_19.rds")) %>% 
    tidysc::aggregate_cells(c(sample, cell_type), slot = "counts")

in_data = in_data %>% 
    mutate(response = fct_recode(severity, 
                                 "Not Critical" = "control",
                                 "Not Critical" = "moderate",
                                 "Critical" = "critical"))

in_data = in_data %>% nest(data_cell_type = -cell_type) %>% 
    mutate(saved = map2(
        data_cell_type, cell_type,
        ~ .x %>%
            mutate(cell_type = .y) %>%
            saveRDS(glue("{out_dir}{.y}_raw.rds"))
    ))