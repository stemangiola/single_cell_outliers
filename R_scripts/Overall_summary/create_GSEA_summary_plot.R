library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ggrepel)
library(ggallin)

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
                  vjust = 1,
                )
            )
    )

in_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/GSEA_data/" # path to data
out_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/"
de_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/"
files = list.files(path = glue("{in_dir}")) # get the file name

df <- files %>% 
    map(~ readRDS(file.path(in_dir, .))) %>% 
    # purrr::reduce(rbind)
    dplyr::bind_rows()

df %>% arrange(overlapped_prop_in_all_pathways)

COVID_19 <- readRDS(glue("{de_dir}COVID_19/data/final_sum_table/all_de.rds"))

COVID_19_sum_counts = COVID_19 %>% 
    # filter(method == !!method) %>% 
    filter(method == "edgeR_quasi_likelihood") %>%
    dplyr::group_by(cell_type, sample) %>% 
    dplyr::summarise(sum_counts = sum(abundance_RNA)) %>% ungroup() 
COVID_19_grand_sum_counts = COVID_19_sum_counts %>% 
    group_by(sample) %>% 
    dplyr::summarise(grand_sum_counts = sum(sum_counts)) %>% ungroup()
COVID_19_sum_counts <- COVID_19_sum_counts %>% left_join(COVID_19_grand_sum_counts, by = "sample")
COVID_19_sum_counts <- COVID_19_sum_counts %>% 
    mutate(proportion_of_all_cell_types = sum_counts/grand_sum_counts) 
COVID_19_average_RNA_proportion <- COVID_19_sum_counts %>% 
    group_by(cell_type) %>% 
    dplyr::summarise(average_RNA_proportion_of_cell_type_among_all_cell_types = mean(proportion_of_all_cell_types))
COVID_19_sum_counts <- COVID_19_sum_counts %>% left_join(COVID_19_average_RNA_proportion, by = "cell_type")

# Extract the column:average_RNA_proportion_of_cell_type_among_all_cell_types
df2 <- COVID_19_sum_counts %>% distinct(cell_type, average_RNA_proportion_of_cell_type_among_all_cell_types)

# # Calculate the proportions of the outliers for each cell type
# COVID_19_summary_panelD <- COVID_19_summary_table %>%
#     pivot_wider(names_from = variable, values_from = Count) %>% 
#     mutate(prop_outlier = genes_with_outliers/de_genes) %>%
#     select(c(cell_type,prop_outlier)) %>%
#     mutate(prop_outlier = replace_na(prop_outlier,0))
# 
# COVID_19_panelD_dataframe <- COVID_19_sum_counts %>%
#     left_join(COVID_19_summary_panelD, by = "cell_type") %>%
#     mutate(dataset_ID = "COVID_19")
df2$cell_type <- gsub(" ", "_", df2$cell_type)

df3 <- df %>% left_join(df2, by = "cell_type")

bar_plot = df3 %>% ggplot(aes(x = cell_type, y = overlapped_prop_in_all_pathways)) +
    geom_point(aes(colour = method, size = average_RNA_proportion_of_cell_type_among_all_cell_types)) +
    # geom_point(aes(colour = method)) +
    custom_theme +
    ggtitle("The plot of the overlapped GSEA pathways proportions including and excluding the outlier transcripts")


