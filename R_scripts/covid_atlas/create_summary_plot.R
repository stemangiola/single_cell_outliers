# Generate 5 summary plots
# aggregate data from all cell types together

library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ggrepel)
library(ggallin)
library(patchwork)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/7e948a2aac3f8ad5989822324395d7d93b51a949/ggplot_theme_multipanel")

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

data_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/data/"
plot_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/covid_atlas/plot/"

# import dataset: merged_summary_table for ppcseq
merged_sum_table <- readRDS(glue("{data_dir}final_sum_table/merged_sum_table.rds"))

# PanelB Summary barplot: x_axis: cell_type | y_axis: proportions of the outliers
merged_sum_table <- merged_sum_table %>% 
    mutate(outlier_prop = num_de_genes_with_outliers/num_de_genes) %>% 
    mutate_at(vars(outlier_prop), ~replace(., is.nan(.), 0))

bar_plot1 = merged_sum_table %>% 
    ggplot(aes(x = reorder(cell_type, -outlier_prop),y = outlier_prop, fill = cell_type)) +
    geom_bar(stat = "identity") +
    scale_y_continuous() +
    custom_theme +
    xlab("cell_type") +
    ylab("proportions of outliers in DE genes") +
    ggtitle("Barplot describing the outlier proportions in DE genes for each cell type") +
    multipanel_theme +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) 

# PanelB stacked summary_barplot: x_axis: cell_type | y_axis: 4 variables
# merged_sum_table <- merged_sum_table %>% select(-outlier_prop)
names(merged_sum_table) <- c("cell_type","total_genes", "de_genes","de_genes_with_outliers", "outliers_in_top_PValue", "outlier_prop")
sum_table = merged_sum_table %>% pivot_longer(cols = c(total_genes, de_genes, de_genes_with_outliers, outliers_in_top_PValue),
                                              names_to = "variable",
                                              values_to = "Count")

sum_table %>% saveRDS(glue("{out_dir}{method}_sum_table_for_bar_plot.rds"))
# draw bar plot

sum_table$variable = factor(sum_table$variable, levels = c("total_genes","de_genes","de_genes_with_outliers","outliers_in_top_PValue"))
dev.new()
bar_plot2 <- sum_table %>% ggplot(aes(x = reorder(cell_type, -outlier_prop), y = Count,fill = variable)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(trans = "log10") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("cell_type") +
    ylab("log(Count)") +
    ggtitle("summary barplot of covid_atlas") 
    # multipanel_theme

library(patchwork)

PanelB <- (bar_plot1 / bar_plot2) + #guide_area() +
    plot_layout(widths = c(0.8,1), guides = 'collect') &
    theme(
        plot.tag = element_text(size = 7),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        legend.key.size = unit(0.1, 'cm')
    )

## Volcano plot
de_files = list.files(path = glue("{data_dir}de_data/")) # get the file name

df <- de_files %>% 
    map(~ readRDS(file.path(glue("{data_dir}de_data/"), .))) %>% 
    # purrr::reduce(rbind)
    dplyr::bind_rows()

# read all ppcseq results
ppcseq_dir = glue("{data_dir}ppcseq_data/")
ppcseq_files = list.files(path = glue("{ppcseq_dir}"))

ppcseq_df <- tibble(filename = ppcseq_files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(ppcseq_dir, .))) # a new data column
    ) 

colnames(ppcseq_df)[1] <- "cell_type"
ppcseq_df$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df$cell_type)
# ppcseq_df %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/all_genes_with_outliers.rds")

ppcseq_df = ppcseq_df %>% unnest(cols = file_contents) %>% unnest(sample_wise_data) %>% filter(deleterious_outliers)

data_for_volanco <- ppcseq_df %>% left_join(df, by = c("cell_type", "transcript", "abundance_RNA", "sample_name"))

data_for_volanco <-  data_for_volanco %>% 
    select(-c(`F___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - disease_severity_standardhealthy`,
              `F___disease_severity_standardsevere - (1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardhealthy)`,
              `PValue___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - disease_severity_standardhealthy`,
              `PValue___disease_severity_standardsevere - (1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardhealthy)`,
              `FDR___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - disease_severity_standardhealthy`,
              `FDR___disease_severity_standardsevere - (1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardhealthy)`,
              `logFC___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - disease_severity_standardhealthy`,
              `logFC___disease_severity_standardsevere - (1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardhealthy)`,
              `logCPM___(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - disease_severity_standardhealthy`,
              `logCPM___disease_severity_standardsevere - (1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardhealthy)`))

# data_for_volanco <- data_for_volanco %>% nest(data = -cell_type) %>%
#     mutate(test_df = map(data, ~.x %>% pivot_transcript(
#         .transcript = transcript
#     ))) %>%
#     select(cell_type,test_df) %>% unnest(test_df)


# data_for_volanco %>% 
#     ggplot(aes(x = `logFC___disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)`,
#                y = `PValue___disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)`,
#                color = cell_type)) +
#     geom_point() +
#     scale_y_continuous(trans = "log10_reverse") +
#     geom_hline(yintercept= 0.05, linetype='dotted', col = 'red') +
#     geom_text(aes(-2, 0.05 ,label = "PValue = 0.05", vjust = -1, colour = "red")) +
#     xlab("logFC") +
#     ylab("PValue") +
#     ggtitle("Volcano plot of the outlier DE genes in the contrast: severe COVID versus (mild and moderate COVID and healthy")+
#     custom_theme 

data_for_volanco <- data_for_volanco %>% 
# Group by contrast. Comparisons both ways.
    pivot_longer(
        cols = contains("___"),
        names_to = c("stats", "contrast"),
        values_to = ".value",
        names_sep="___"
    ) %>%

    # Markers selection within each pair of contrast
    nest(stat_df = -contrast) %>%

    # Reshape inside each contrast
    mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
    unnest(stat_df)

data_for_volanco <- data_for_volanco %>% nest(data = -contrast) %>% 
    mutate(data = map(data, ~ .x %>% filter(PValue < 0.05))) %>% unnest(data)

contrast_names <- c(
    "(1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardsevere) - disease_severity_standardhealthy" = "covid vs healthy",
    "(1/2 * disease_severity_standardmoderate + 1/2 * disease_severity_standardsevere) - ( 1/2 * disease_severity_standardhealthy + 1/2 * disease_severity_standardmild)" = "(moderate & severe) vs (healthy & mild)",
    "disease_severity_standardsevere - (1/3 * disease_severity_standardmild + 1/3 * disease_severity_standardmoderate + 1/3 * disease_severity_standardhealthy)" = "severe vs (heathy & mild & moderate)"
)

data_for_volanco %>% 
    ggplot(aes(x = logFC,
               y = PValue,
               color = cell_type)) +
    geom_point() +
    scale_y_continuous(trans = "log10_reverse") +
    geom_hline(yintercept= 0.05, linetype='dotted', col = 'red') +
    geom_text(aes(-2, 0.05 ,label = "PValue = 0.05", vjust = -1, colour = "red")) +
    xlab("logFC") +
    ylab("PValue") +
    facet_wrap(~ contrast, labeller = as_labeller(contrast_names)) +
    ggtitle("Volcano plot of the outlier DE genes in 3 contrasts")+
    multipanel_theme +
    custom_theme 



    
# PanelD: Overlapped histogram Summary histogram（Histogram_fold_change_plot）
# read all ppcseq results
ppcseq_dir = glue("{data_dir}ppcseq_data/")
ppcseq_files = list.files(path = glue("{ppcseq_dir}"))

ppcseq_df <- tibble(filename = ppcseq_files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(ppcseq_dir, .))) # a new data column
    ) 

colnames(ppcseq_df)[1] <- "cell_type"
ppcseq_df_FC = ppcseq_df %>% unnest(cols = file_contents) %>% 
    select(cell_type, transcript, sample_wise_data) %>% 
    unnest(cols = sample_wise_data) 

# change cell_type from ....rds to ""
ppcseq_df_FC$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df_FC$cell_type)
# ppcseq_df$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df$cell_type)

# generate the df for calculating the FC before and after outlier exclusion
ppcseq_df_FC = ppcseq_df_FC %>% 
    filter(deleterious_outliers) %>% 
    distinct(cell_type, transcript, slope_after_outlier_filtering, slope_before_outlier_filtering) %>% 
    mutate(FC = slope_after_outlier_filtering/slope_before_outlier_filtering) 
ppcseq_df_FC %>% 
    saveRDS(glue("{data_dir}final_sum_table/FC_df.rds"))

FC_plot <- ppcseq_df_FC %>% ggplot(aes(x = FC, color =cell_type)) +
    # geom_histogram(bins = 30) +
    geom_density() +
    scale_x_continuous() +
    xlim(-2, 3) +
    # facet_wrap(~ cell_type, scales = "free") + 
    # ggtitle("Histogram of fold changes of differential transcripts with and without outliers for method: deseq2") +
    ggtitle(glue("Ratio of fold changes of differential transcripts after excluding the outliers")) +
    xlab("ratio of fold change after outlier elimination") + 
    custom_theme


## Density
dev.new()
density_plot <- df %>% 
    ggplot(aes(x = abundance_RNA_scaled + 1, color = sample)) +
    geom_density() +
    facet_wrap(~ cell_type) +
    scale_x_log10() +
    ggtitle("Density plot of 31 cell types of covid_atlas") +
    theme(legend.position="none")

density_plot %>% ggsave(filename = "Density_plot.pdf", 
                        path = glue("{plot_dir}"),
                        width = 14, height = 10, units = "in", dpi = 300)

dev.off()
#####------------------------------------------------------------------------------------------------------------
# # PCA plot
# df_for_PCA = df %>% nest(data = -cell_type) %>% 
#     mutate(test_df = map(data, ~.x %>% pivot_sample())) %>% 
#     select(cell_type,test_df) %>% unnest(test_df)
# dev.new()
# PCA_plot <- df_for_PCA %>% 
#     # ggplot(aes(x = PC1, y = PC2, colour = log(size))) +
#     ggplot(aes(x = PC1, y = PC2, colour = disease_severity_standard)) +
#     # scale_fill_distiller(palette = "Spectral") +
#     geom_point() +
#     facet_wrap(~ cell_type) +
#     ggtitle(glue("PCA plot for 31 cell types covid_atlas")) +
#     custom_theme
# 
# # Save PCA_plot as png & pdf
# PCA_plot %>% ggsave(filename = glue("{method}_PCA_plot.pdf"), 
#                     path = glue("{plot_dir}"),
#                     width = 14, height = 10, units = "in", dpi = 300)
# 
# PCA_plot %>% ggsave(filename =  glue("{method}_PCA_plot.png"), 
#                     path = glue("{plot_dir}"),
#                     width = 14, height = 10, units = "in", dpi = 300)
# dev.off()


