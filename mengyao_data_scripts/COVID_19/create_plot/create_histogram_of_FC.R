# test FC with and without outliers in DE genes
library(tidyverse)
library(glue)
library(ggallin)
library(ggplot2)

# df =  ppcseq.rds
ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/ppcseq_data/"
ppcseq_files = list.files(path = glue("{ppcseq_dir}"))
cell_names = tools::file_path_sans_ext(ppcseq_files) %>% basename() %>% str_replace_all("_ppcseq","")

ppcseq_df <- tibble(filename = ppcseq_files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(ppcseq_dir, .))) # a new data column
    ) 

colnames(ppcseq_df)[1] <- "cell_type"
# change cell_type from ....rds to ""
ppcseq_df$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df$cell_type)

# generate the df for calculating the FC before and after outlier exclusion
ppcseq_df = ppcseq_df %>% unnest(cols = file_contents) %>% 
    select(cell_type, transcript, sample_wise_data) %>% 
    unnest(cols = sample_wise_data) 

ppcseq_df = ppcseq_df %>% 
    distinct(cell_type, transcript, slope_after_outlier_filtering, slope_before_outlier_filtering) %>% 
    mutate(FC = slope_after_outlier_filtering/slope_before_outlier_filtering) 

ppcseq_df %>% ggplot(aes(x = FC)) +
    geom_histogram(bins = 30) +
    xlim(-10, 10) +
    scale_x_continuous(trans = pseudolog10_trans, guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ cell_type, scales = "free") + 
    ggtitle("Histogram of fold changes of differential transcripts with and without outliers") +
    custom_theme 
