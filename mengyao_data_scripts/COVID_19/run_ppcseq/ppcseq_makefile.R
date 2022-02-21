# ppcseq_makefile.R
# Create makefile.R
library(glue)
library(tidyverse)

script_dir = "/stornext/HPCScratch/home/ma.m/mengyao_data_scripts/COVID_19/run_ppcseq/"
input_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/de_data/"
out_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/ppcseq_data/"
# error_file = "> /stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/ppcseq_data/err.stderr  2>&1"

files = list.files(path = glue("{input_dir}"))
cell_name = tools::file_path_sans_ext(files) %>% basename() %>% str_replace_all("_DE","")

command_df = tibble(out_dir = out_dir, 
                    cell_name = cell_name,
                    Rscript = glue("Rscript {script_dir}run_ppcseq.R"),
                    input_dir = input_dir,
                    files = files) %>% 
  unite("output_file", c(out_dir, cell_name), sep = "") %>% 
  mutate(output_file = glue("{output_file}_ppcseq.rds")) %>% 
  unite("args", c(input_dir, files), sep = "") %>% 
  mutate(command = sprintf("%s:\n\t%s %s %s",
                           output_file, Rscript, args, output_file))

# pull command and add SLURM requirements
command_df %>% pull(command) %>% 
  prepend("CATEGORY=run_ppcseq\nMEMORY=30024\nCORES=4") %>% 
  write_lines(glue("{script_dir}run_ppcseq.makeflow"))