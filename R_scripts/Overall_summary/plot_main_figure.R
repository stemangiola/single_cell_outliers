library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ppcseq)
library(plyr)
# devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
friendly_cols <- dittoSeq::dittoColors()
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
                text = element_text(size = 16),
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
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/7e948a2aac3f8ad5989822324395d7d93b51a949/ggplot_theme_multipanel")

in_dir <- "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/"
# method = c("deseq2", "edgeR_quasi_likelihood","edger_robust_likelihood_ratio")
# files = list.files(path = glue("{in_dir}COVID_19/data/final_sum_table/"), pattern = "_final_sum_table.rds") # get the file name
method = "edgeR_quasi_likelihood"
# multipanel_theme
# Summary Panel A 
# Use final summary table 
COVID_19_summary_table <- list.files(path = glue("{in_dir}COVID_19/data/final_sum_table/"), 
                                      pattern = "_final_sum_table.rds") %>% 
    map(~ readRDS(file.path(glue("{in_dir}COVID_19/data/final_sum_table/"), .))) %>% 
    dplyr::bind_rows()
COVID_19_summary_panelA <- COVID_19_summary_table %>% 
    pivot_wider(names_from = variable, values_from = Count) %>% 
    select(-c(transcript, outliers_in_top_PValue)) %>% 
    dplyr::group_by(method) %>%
    # filter(method == method) %>%
    dplyr::summarise(genes_with_outliers = sum(genes_with_outliers), de_genes = sum(de_genes)) %>% 
    mutate(prop_outliers = genes_with_outliers/de_genes) %>% mutate(name = "COVID_19")


GSE120575_summary_table <- list.files(path = glue("{in_dir}GSE120575_melanoma/data/final_sum_table/"), 
                                      pattern = "_sum_table_for_bar_plot.rds") %>% 
    map(~ readRDS(file.path(glue("{in_dir}GSE120575_melanoma/data/final_sum_table/"), .))) %>% 
    dplyr::bind_rows()
GSE120575_summary_panelA <- GSE120575_summary_table %>% 
    pivot_wider(names_from = variable, values_from = Count) %>% 
    select(-c(transcript, outliers_in_top_PValue)) %>% 
    dplyr::group_by(method) %>%
    # filter(method == method) %>% 
    dplyr::summarise(genes_with_outliers = sum(genes_with_outliers), de_genes = sum(de_genes)) %>% 
    mutate(prop_outliers = genes_with_outliers/de_genes) %>% mutate(name = "GSE120575")


GSE129788_summary_table <- list.files(path = glue("{in_dir}GSE129788_aging_mouse_brain/data/final_sum_table/"), 
                                      pattern = "_sum_table_for_bar_plot.rds") %>% 
    map(~ readRDS(file.path(glue("{in_dir}GSE129788_aging_mouse_brain/data/final_sum_table/"), .))) %>% 
    dplyr::bind_rows()
GSE129788_summary_panelA <- GSE129788_summary_table %>% 
    pivot_wider(names_from = variable, values_from = Count) %>% 
    select(-c(transcript, outliers_in_top_PValue)) %>% 
    dplyr::group_by(method) %>%
    # filter(method == method) %>%
    dplyr::summarise(genes_with_outliers = sum(genes_with_outliers), de_genes = sum(de_genes)) %>% 
    mutate(prop_outliers = genes_with_outliers/de_genes) %>% mutate(name = "GSE129788")

SCP1288_summary_table <- list.files(path = glue("{in_dir}SCP1288_renal_cell_carcinoma/data/final_sum_table/"), 
                                    pattern = "_sum_table_for_bar_plot.rds") %>% 
    map(~ readRDS(file.path(glue("{in_dir}SCP1288_renal_cell_carcinoma/data/final_sum_table/"), .))) %>% 
    dplyr::bind_rows()
SCP1288_summary_panelA <- SCP1288_summary_table %>% 
    pivot_wider(names_from = variable, values_from = Count) %>% 
    select(-c(transcript, outliers_in_top_PValue)) %>%
    dplyr::group_by(method) %>%
    # filter(method == method) %>% 
    dplyr::summarise(genes_with_outliers = sum(genes_with_outliers), de_genes = sum(de_genes)) %>% 
    mutate(prop_outliers = genes_with_outliers/de_genes) %>% mutate(name = "SCP1288")

PanelA_df <- COVID_19_summary_panelA %>% rbind(GSE120575_summary_panelA, GSE129788_summary_panelA, SCP1288_summary_panelA) 


# Barplot
PanelA_plot <- PanelA_df %>% 
    # filter(method == "edgeR_quasi_likelihood") %>% 
    ggplot(aes(x=name, y=prop_outliers, fill = method)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    ylab("Fraction of significant transcripts including outliers") +
    # ylab(glue("Fraction of Significant transcripts including outliers with method")) +
    xlab("Data set") +
    multipanel_theme + 
    custom_theme

# Panel B ppcseq plot -----------------------------------------------------------------------------------------------------
# plot_credible_intervals() %>% pull(plot) %>% .[1]

# input datasets -- COVID_19_summary_table
# select cell_type pDC
COVID_19_pDC <- readRDS(glue("{in_dir}COVID_19/data/ppcseq_data/{method}_pDC_ppcseq.rds"))
# select gene IGJ
PanelB_1 <- COVID_19_pDC %>% plot_credible_intervals() %>% filter(transcript == "IRF2BP2") %>% pull(plot) 
PanelB_1 <- PanelB_1[[1]]

# select cell_type, Lymphocytes_exhausted_cell_cycle
GSE120575_Exhausted_CD8_T_cells <- 
    readRDS(glue("{in_dir}GSE120575_melanoma/data/ppcseq_data/{method}_Exhausted_CD8+_T_cells_ppcseq.rds"))
PanelB_2 <- GSE120575_Exhausted_CD8_T_cells %>% plot_credible_intervals() %>% 
    filter(transcript == "ATP1B3") %>% pull(plot)
# PanelB_2 <- PanelB_2[[1]] + multipanel_theme
PanelB_2 <- PanelB_2[[1]] 

GSE129788_Ependymocytes <- 
    readRDS(glue("{in_dir}GSE129788_aging_mouse_brain/data/ppcseq_data/{method}_Ependymocytes_ppcseq.rds"))
PanelB_3 <- GSE129788_Ependymocytes %>% plot_credible_intervals() %>% 
    filter(transcript == "H2-D1") %>% pull(plot) 
PanelB_3 <- PanelB_3[[1]]


SCP1288_Misc_Undetermined <- readRDS(glue("{in_dir}SCP1288_renal_cell_carcinoma/data/ppcseq_data/{method}_Misc_Undetermined_ppcseq.rds"))
PanelB_4 <- SCP1288_Misc_Undetermined %>% plot_credible_intervals() %>%
    filter(transcript == "TAGLN") %>% pull(plot) 
PanelB_4 <- PanelB_4[[1]]


# Panel C:density plot for the differences in fold change caused by outlier elimination -----------------------------------------------------------------------------------------------------
# input datasets
COVID19_FC <- readRDS(glue("{in_dir}COVID_19/data/final_sum_table/{method}_FC_df.rds"))
COVID19_FC <- COVID19_FC %>% mutate(dataset_ID = "COVID_19")

GSE120575_FC <- readRDS(glue("{in_dir}GSE120575_melanoma/data/final_sum_table/{method}_FC_df.rds"))
GSE120575_FC <- GSE120575_FC %>% mutate(dataset_ID = "GSE120575")

GSE129788_FC <- readRDS(glue("{in_dir}GSE129788_aging_mouse_brain/data/final_sum_table/{method}_FC_df.rds"))
GSE129788_FC <- GSE129788_FC %>% mutate(dataset_ID = "GSE129788")

SCP1288_FC <- readRDS(glue("{in_dir}SCP1288_renal_cell_carcinoma/data/final_sum_table/{method}_FC_df.rds"))
SCP1288_FC  <- SCP1288_FC %>% mutate(dataset_ID = "SCP1288")

# ratio of fold change
# difference of fold change

# merge 4 datasets together
PanelC_df <- COVID19_FC  %>% rbind(GSE120575_FC, GSE129788_FC, SCP1288_FC) 


# draw Fold change density plot
PanelC_plot = PanelC_df %>% ggplot(aes(x = FC, color = cell_type)) +
    geom_density() +
    geom_vline(xintercept = 1, linetype="dotted", 
               color = "grey", size=1.5) +
    xlim(-10, 10) +
    facet_wrap(~ dataset_ID, 
               scales = "free") +
    xlab("Ratio of fold change") +
    ylab("Density") +
    ggtitle(glue("Ratio of fold change after outlier elimination with method: {method}")) +
    multipanel_theme +
    custom_theme

# Panel D scatter plot of average (across samples) of the  total sample RNA vs proportion of outliers in significant genes--------------------------------------------
# Import dataset COVID_19 DE_results
COVID_19 <- readRDS(glue("{in_dir}COVID_19/data/final_sum_table/all_de.rds"))

COVID_19_sum_counts = COVID_19 %>% 
    filter(method == !!method) %>% 
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

# # Calculate the proportions of the outliers for each cell type
COVID_19_summary_panelD <- COVID_19_summary_table %>%
    pivot_wider(names_from = variable, values_from = Count) %>% 
    mutate(prop_outlier = genes_with_outliers/de_genes) %>%
    select(c(cell_type,prop_outlier)) %>%
    mutate(prop_outlier = replace_na(prop_outlier,0))

COVID_19_panelD_dataframe <- COVID_19_sum_counts %>%
    left_join(COVID_19_summary_panelD, by = "cell_type") %>%
    mutate(dataset_ID = "COVID_19")

# GSE120575
GSE120575 <- readRDS(glue("{in_dir}GSE120575_melanoma/data/final_sum_table/all_de.rds"))
GSE120575_sum_counts = GSE120575 %>%
    filter(method == !!method) %>% 
    dplyr::group_by(cell_type, sample) %>% 
    dplyr::summarise(sum_counts = sum(abundance_RNA)) %>% ungroup() 

GSE120575_grand_sum_counts = GSE120575_sum_counts %>% 
    group_by(sample) %>% 
    dplyr::summarise(grand_sum_counts = sum(sum_counts)) %>% ungroup()
GSE120575_sum_counts <- GSE120575_sum_counts %>% left_join(GSE120575_grand_sum_counts, by = "sample")
GSE120575_sum_counts <- GSE120575_sum_counts %>% 
    mutate(proportion_of_all_cell_types = sum_counts/grand_sum_counts) 
GSE120575_average_RNA_proportion <- GSE120575_sum_counts %>% 
    group_by(cell_type) %>% 
    dplyr::summarise(average_RNA_proportion_of_cell_type_among_all_cell_types = mean(proportion_of_all_cell_types))
GSE120575_sum_counts <- GSE120575_sum_counts %>% left_join(GSE120575_average_RNA_proportion, by = "cell_type")

# # Calculate the proportions of the outliers for each cell type
GSE120575_summary_panelD <- GSE120575_summary_table %>%
    pivot_wider(names_from = variable, values_from = Count) %>%
    mutate(prop_outlier = genes_with_outliers/de_genes) %>%
    select(c(cell_type,prop_outlier)) %>%
    mutate(prop_outlier = replace_na(prop_outlier,0))

GSE120575_panelD_dataframe <- GSE120575_sum_counts %>%
    left_join(GSE120575_summary_panelD, by = "cell_type") %>%
    mutate(dataset_ID = "GSE120575")

   
# GSE129788
GSE129788 <- readRDS(glue("{in_dir}GSE129788_aging_mouse_brain/data/final_sum_table/all_de.rds"))

GSE129788_sum_counts = GSE129788 %>% 
    filter(method == !!method) %>% 
    dplyr::group_by(cell_type, sample) %>% 
    dplyr::summarise(sum_counts = sum(abundance_RNA)) %>% ungroup() 

GSE129788_grand_sum_counts = GSE129788_sum_counts %>% 
    group_by(sample) %>% 
    dplyr::summarise(grand_sum_counts = sum(sum_counts)) %>% ungroup()
GSE129788_sum_counts <- GSE129788_sum_counts %>% left_join(GSE129788_grand_sum_counts, by = "sample")
GSE129788_sum_counts <- GSE129788_sum_counts %>% 
    mutate(proportion_of_all_cell_types = sum_counts/grand_sum_counts) 
GSE129788_average_RNA_proportion <- GSE129788_sum_counts %>% 
    group_by(cell_type) %>% 
    dplyr::summarise(average_RNA_proportion_of_cell_type_among_all_cell_types = mean(proportion_of_all_cell_types))
GSE129788_sum_counts <- GSE129788_sum_counts %>% left_join(GSE129788_average_RNA_proportion, by = "cell_type")


# # Calculate the proportions of the outliers for each cell type
GSE129788_summary_panelD <- GSE129788_summary_table %>%
    pivot_wider(names_from = variable, values_from = Count) %>%
    mutate(prop_outlier = genes_with_outliers/de_genes) %>%
    select(c(cell_type,prop_outlier)) %>%
    mutate(prop_outlier = replace_na(prop_outlier,0))

GSE129788_panelD_dataframe <- GSE129788_sum_counts %>%
    left_join(GSE129788_summary_panelD, by = "cell_type") %>%
    mutate(dataset_ID = "GSE129788")


# SCP1288
SCP1288 <- readRDS(glue("{in_dir}SCP1288_renal_cell_carcinoma/data/final_sum_table/all_de.rds"))
colnames(SCP1288)
colnames(SCP1288)[21] <- "abundance_RNA"
colnames(SCP1288)[28] <- "abundance_RNA_scaled"

SCP1288_sum_counts = SCP1288 %>% 
    filter(method == !!method) %>%
    dplyr::group_by(cell_type, sample) %>%
    dplyr::summarise(sum_counts = sum(abundance_RNA)) %>% ungroup()

SCP1288_grand_sum_counts = SCP1288_sum_counts %>% 
    group_by(sample) %>% 
    dplyr::summarise(grand_sum_counts = sum(sum_counts)) %>% ungroup()
SCP1288_sum_counts <- SCP1288_sum_counts %>% left_join(SCP1288_grand_sum_counts, by = "sample")
SCP1288_sum_counts <- SCP1288_sum_counts %>% 
    mutate(proportion_of_all_cell_types = sum_counts/grand_sum_counts) 
SCP1288_average_RNA_proportion <- SCP1288_sum_counts %>% 
    group_by(cell_type) %>% 
    dplyr::summarise(average_RNA_proportion_of_cell_type_among_all_cell_types = mean(proportion_of_all_cell_types))
SCP1288_sum_counts <- SCP1288_sum_counts %>% left_join(SCP1288_average_RNA_proportion, by = "cell_type")


# # Calculate the proportions of the outliers for each cell type
SCP1288_summary_panelD <- SCP1288_summary_table %>%
    pivot_wider(names_from = variable, values_from = Count) %>%
    mutate(prop_outlier = genes_with_outliers/de_genes) %>%
    select(c(cell_type,prop_outlier)) %>%
    mutate(prop_outlier = replace_na(prop_outlier,0))

SCP1288_panelD_dataframe <- SCP1288_sum_counts %>%
    left_join(SCP1288_summary_panelD, by = "cell_type") %>%
    mutate(dataset_ID = "SCP1288")


# COVID_19 %>% 
#     nest(data = -c(sample,cell_type, dataset_id)) %>% 
#     mutate(sum_counts = map_dbl(data, ~.x$abundance_RNA %>% sum())) %>% 
#     with_groups(sample, ~ mutate(.x, grand_sum_counts = sum(sum_counts))) %>% 
#     mutate(proportion_of_all_cell_types = sum_counts/grand_sum_counts) %>% 
#     with_groups(cell_type,  ~ summarise(.x, average_RNA_proportion_of_cell_type_among_all_cell_types = mean(proportion_of_all_cell_types)))

PanelD_df <- COVID_19_panelD_dataframe %>% rbind(GSE120575_panelD_dataframe, GSE129788_panelD_dataframe, SCP1288_panelD_dataframe) 

PanelD_plot <- PanelD_df %>% ggplot(aes(x = average_RNA_proportion_of_cell_type_among_all_cell_types, 
                                        y = prop_outlier, size = grand_sum_counts)) +
    geom_point() +
    facet_wrap(~ dataset_ID, 
               scales = "free") +
    xlab("Average RNA proportions among all cell types") +
    ylab("Outlier proportions in significant genes") +
    multipanel_theme + custom_theme


# Panel E: Venn diagram that counts how many genes including outliers overlap across all three methods.

COVID_19_df_for_venn <- list.files(path = glue("{in_dir}COVID_19/data/ppcseq_data/")) %>% 
    map(~ readRDS(file.path(glue("{in_dir}COVID_19/data/ppcseq_data/"), .))) %>% 
    dplyr::bind_rows() 
COVID_19_df_for_venn <- COVID_19_df_for_venn %>% 
    select(c(transcript, method, tot_deleterious_outliers)) %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    mutate(dataset_ID = "COVID_19")

GSE120575_df_for_venn <- list.files(path = glue("{in_dir}GSE120575_melanoma/data/ppcseq_data/")) %>% 
    map(~ readRDS(file.path(glue("{in_dir}GSE120575_melanoma/data/ppcseq_data/"), .))) %>% 
    dplyr::bind_rows() 
    
GSE120575_df_for_venn <- GSE120575_df_for_venn %>% 
    select(c(transcript, method, tot_deleterious_outliers)) %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    mutate(dataset_ID = "GSE120575")

GSE129788_df_for_venn <- list.files(path = glue("{in_dir}GSE129788_aging_mouse_brain/data/ppcseq_data/")) %>% 
    map(~ readRDS(file.path(glue("{in_dir}GSE129788_aging_mouse_brain/data/ppcseq_data/"), .))) %>% 
    dplyr::bind_rows() 
GSE129788_df_for_venn <- GSE129788_df_for_venn %>% 
    select(c(transcript, method, tot_deleterious_outliers)) %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    mutate(dataset_ID = "GSE129788")

SCP1288_df_for_venn <- list.files(path = glue("{in_dir}SCP1288_renal_cell_carcinoma/data/ppcseq_data/")) %>% 
    map(~ readRDS(file.path(glue("{in_dir}SCP1288_renal_cell_carcinoma/data/ppcseq_data/"), .))) %>% 
    dplyr::bind_rows() 

SCP1288_df_for_venn <- SCP1288_df_for_venn %>% 
    select(c(transcript, method, tot_deleterious_outliers)) %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    mutate(dataset_ID = "SCP1288")

PanelE_df <- COVID_19_df_for_venn %>% rbind(GSE120575_df_for_venn, GSE129788_df_for_venn, SCP1288_df_for_venn)
PanelE_df <- PanelE_df %>% 
    unite("transcript_name", c(dataset_ID, transcript), sep = "_", remove = FALSE)

deseq2_list <- PanelE_df %>% filter(method == "deseq2") %>% pull(transcript_name) 
edgeR_quasi_likelihood_list <- PanelE_df %>% filter(method == "edgeR_quasi_likelihood") %>% pull(transcript_name) 
edger_robust_likelihood_ratio_list <- PanelE_df %>% filter(method == "edger_robust_likelihood_ratio") %>% pull(transcript_name)

venn_diagram_list <- list(deseq2 = deseq2_list, 
                          edgeR = edgeR_quasi_likelihood_list, 
                          edgeR_robust = edger_robust_likelihood_ratio_list)
PanelE_plot <- ggvenn(
    venn_diagram_list, 
    fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
    stroke_size = 0.4, set_name_size = 3,
    text_size = 3
)


# Generate the main figure 
library(patchwork) 
# 20 * 15 inch
p = (
    (PanelA_plot - (PanelB_1 + PanelB_2 + PanelB_3 + PanelB_4 + plot_layout(heights = c(1, 1, 1, 1), nrow = 2)) +
         plot_layout(widths = c(0.8,2), guides = 'keep')) /
        (PanelC_plot + 
             (PanelD_plot + plot_spacer() + PanelE_plot + plot_layout(widths = c(2, 0.2,1.5)))) 
) + plot_layout(heights = c(2.5,1.2),nrow = 2) + 
    plot_annotation(tag_levels = c('A')) &
    theme(
        plot.tag = element_text(size = 8),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "pt"),
        legend.key.size = unit(0.1, 'cm')
    )

# p =
#     (
#         (
#             plot_UMAP_samples + plot_UMAP_all + plot_generic_markers + plot_bar_generic +
#                 plot_layout(height = c(1, 1, 1, 1), nrow = 1)
#         )  /
#             (
#                 plot_UMAP_lymphoid + plot_lymphoid_markers + plot_spacer() + plot_layout(widths = c(1, 2.5, 0.8), nrow = 1)
#             ) /
#             plot_bar_lymphoid /
#             
#             (
#                 plot_UMAP_myeloid + plot_myeloid_markers + plot_spacer() + plot_layout(widths = c(1, 2.5, 0.8), nrow = 1)
#             ) /
#             plot_bar_myeloid
#     ) +
#     plot_layout(height  = c(1, 2, 1, 2, 1), nrow = 5)  +
#     plot_annotation(tag_levels = c('A')) &
    # theme(
    #     plot.tag = element_text(size = 8),
    #     plot.margin = margin(0, 0, 0, 0, "pt"),
    #     legend.key.size = unit(0.2, 'cm')
    # )
















