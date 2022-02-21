library(tidyverse)
library(glue)

script_dir = "/stornext/HPCScratch/home/ma.m/mengyao_data_scripts/COVID_19/create_plot/"
ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/ppcseq_data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/"

ppcseq_files = list.files(path = glue("{ppcseq_dir}"))
cell_name = tools::file_path_sans_ext(ppcseq_files) %>% basename() %>% str_replace_all("_ppcseq","")

command_df = tibble(out_dir = out_dir, 
                    cell_name = cell_name,
                    Rscript = glue("Rscript {script_dir}transcript_list.R"),
                    ppcseq_dir = ppcseq_dir,
                    ppcseq_files = ppcseq_files) %>% 
  unite("output_file", c(out_dir, cell_name), sep = "") %>% 
  mutate(output_file = glue("{output_file}_transcript_list.rds")) %>% 
  unite("arg1", c(ppcseq_dir, ppcseq_files), sep = "") %>% 
  mutate(command = sprintf("%s:\n\t%s %s %s %s",
                           output_file, Rscript, arg1, cell_name, output_file))
# pull command and add SLURM requirements
command_df %>% pull(command) %>% 
  prepend("CATEGORY=create_transcript_list\nMEMORY=30024\nCORES=4") %>% 
  write_lines(glue("{script_dir}create_transcript.makeflow"))