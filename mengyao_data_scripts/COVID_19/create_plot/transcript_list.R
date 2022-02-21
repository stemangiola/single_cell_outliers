library(tidyverse)
library(glue)

args = commandArgs(trailingOnly = TRUE)
filename = args[1]
cell_name = args[2]
out_file = args[3]

input_file = readRDS(file = filename)
cell_name = "Squamous"
if (nrow(input_file != 0)){
  input_file %>% 
    select(transcript, tot_deleterious_outliers) %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    mutate(cell_type = glue("{cell_name}")) %>% 
    nest(data = -cell_type) %>% saveRDS(file = out_file)
} else{
  df = tibble(cell_type = glue("{cell_name}"),
              transcript = character(),
              tot_deleterious_outliers = integer()
              ) %>% nest(data = -cell_type) 
}



  
  