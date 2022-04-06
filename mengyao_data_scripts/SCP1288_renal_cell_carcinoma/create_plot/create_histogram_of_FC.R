# test FC with and without outliers in DE genes
library(tidyverse)
library(glue)
library(ggallin)
library(ggplot2)
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
    list(
        scale_fill_manual(values = friendly_cols),
        scale_color_manual(values = friendly_cols),
        theme_bw() +
            theme(
                panel.border = element_blank(),
                axis.line = element_line(),
                panel.grid.major = element_line(size = 0.2),
                panel.grid.minor = element_line(size = 0.1),
                text = element_text(size = 12),
                legend.position = "bottom",
                strip.background = element_blank(),
                axis.title.x = element_text(margin = margin(
                    t = 10,
                    r = 10,
                    b = 10,
                    l = 10
                )),
                axis.title.y = element_text(margin = margin(
                    t = 10,
                    r = 10,
                    b = 10,
                    l = 10
                )),
                axis.text.x = element_text(
                    angle = 30,
                    hjust = 1,
                    vjust = 1
                )
            )
    )


# df =  ppcseq.rds
ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/data/ppcseq_data/"

ppcseq_files = list.files(path = glue("{ppcseq_dir}"))

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

dev.new()
histogram_plot <- ppcseq_df %>% ggplot(aes(x = FC)) +
    geom_histogram(bins = 30) +
    xlim(-10, 10) +
    scale_x_continuous(trans = pseudolog10_trans, guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ cell_type, scales = "free") + 
    ggtitle("Histogram of fold changes of differential transcripts with and without outliers for dataset (SCP1288)") +
    custom_theme

histogram_plot %>% ggsave(filename = "Histogram_of_FC.pdf", 
                          path = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/plot/",
                          width = 14, height = 10, units = "in", dpi = 300)

histogram_plot %>% ggsave(filename = "Histogram_of_FC.png", 
                          path = "/stornext/HPCScratch/home/ma.m/single_cell_database/SCP1288_renal_cell_carcinoma/plot/",
                          width = 14, height = 10, units = "in", dpi = 300)

dev.off()
