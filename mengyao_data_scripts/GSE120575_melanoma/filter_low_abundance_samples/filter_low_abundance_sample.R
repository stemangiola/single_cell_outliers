# filter out the low abundance samples from the tidy data
library(tidyverse)
library(glue)

args = commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]

# commands to filter out low abundance samples into cell_type_low_abundance_sample.rds
dataframe = readRDS(in_file)

dataframe <- dataframe %>% 
    nest(data = -sample) %>% 
    mutate(size = map_dbl(data, ~ .x$abundance_RNA %>% sum())) %>%
    mutate(ratio_from_largest = size/median(size)) %>%
    filter(ratio_from_largest > 0.1) %>% 
    unnest (cols = data) %>% 
    saveRDS(glue("{out_file}"))