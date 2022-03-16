# aggregate data from all cell types together
library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ggrepel)
# devtools::install_github("slowkow/ggrepel")


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


in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/de_data/" # path to data
files = list.files(path = glue("{in_dir}")) # get the file name

df <- files %>% 
    map(~ readRDS(file.path(in_dir, .))) %>% 
    reduce(rbind)
# df %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/all_de.rds")

# generate faced density plot
dev.new()
density_plot <- df %>% 
    ggplot(aes(x = abundance_RNA_scaled + 1, color = sample)) +
    geom_density() +
    facet_wrap(~ cell_type) +
    scale_x_log10() +
    ggtitle("Density plot for 11 cell types of dataset (GSE120575)") +
    custom_theme 

density_plot %>% ggsave(filename = "Density_plot.pdf", 
                    path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                    width = 14, height = 10, units = "in", dpi = 300)

density_plot %>% ggsave(filename = "Density_plot.png", 
                    path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                    width = 14, height = 10, units = "in", dpi = 300)
dev.off()


# PCA plot
df_for_PCA = df %>% nest(data = -cell_type) %>% 
    mutate(test_df = map(data, ~.x %>% pivot_sample())) %>% 
    select(cell_type,test_df) %>% unnest(test_df)

dev.new()
PCA_plot <- df_for_PCA %>% 
    ggplot(aes(x = PC1, y = PC2, colour = condition)) +
    geom_point() +
    # geom_text_repel(aes(label = sample), show.legend = FALSE) +
    facet_wrap(~ cell_type) +
    ggtitle("PCA plot for 11 cell types for dataset (GSE120575)") +
    custom_theme

PCA_plot %>% ggsave(filename = "PCA_plot.pdf", 
                        path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                        width = 14, height = 10, units = "in", dpi = 300)

PCA_plot %>% ggsave(filename = "PCA_plot.png", 
                        path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                        width = 14, height = 10, units = "in", dpi = 300)
dev.off()

#####------------------------------------------------------------------------------------------------------------
# df1 %>% 
#     #nest_on_cell_type
#     #produce_plot
#     ggplot(aes(x = PC1, y = PC2, colour = response)) +
#     geom_point() +
#     geom_text(aes(label = sample)) +
#     facet_wrap(~ cell_type) +
#     ggtitle("PCA plot of 24 cell types") +
#     custom_theme
# # Visualise cell type by cell type separately, zoom out, decrease the size of text
# # see if all outliers appear to have the same or overlapping sample labels

###-----------------------------------------------------------------------------------------------------------------

# bar plot
df <- df %>% 
    mutate(abundance_RNA = abundance_RNA %>% as.integer) %>% 
    mutate(is_significant = FDR < 0.05) 

new_df = df %>% nest(data = -cell_type) %>% 
    mutate(tot_transcript = map_dbl(data, ~ .x$transcript %>% n_distinct())) %>% 
    select(-data) 

# read final summary table
sum_table <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/final_sum_table/sum_table.rds")


sum_table = new_df %>% left_join(sum_table) 
sum_table[is.na(sum_table)] <- 0
names(sum_table) <- c("cell_type","transcript", "genes_with_outliers","de_genes","outliers_in_top_PValue")
sum_table = sum_table %>% select(cell_type, transcript, de_genes, genes_with_outliers, outliers_in_top_PValue)
sum_table <- sum_table %>% pivot_longer(cols = c(transcript, de_genes, genes_with_outliers, outliers_in_top_PValue),
                                        names_to = "variable",
                                        values_to = "Count")
sum_table %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/final_sum_table/sum_table_for_bar_plot.rds")
# draw bar plot

# change the levels of 4 variables 
sum_table$variable = factor(sum_table$variable, levels = c("transcript","de_genes","genes_with_outliers","outliers_in_top_PValue"))
dev.new()
bar_plot = sum_table %>% ggplot(aes(x = variable, y = Count,fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Count)) +
    scale_y_continuous() +
    custom_theme +
    facet_wrap(~ cell_type) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("summary table of dataset (GSE120575)")

bar_plot %>% ggsave(filename = "summary_bar_plot.pdf", 
                    path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                    width = 14, height = 10, units = "in", dpi = 300)

bar_plot %>% ggsave(filename = "summary_bar_plot.png", 
                    path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                    width = 14, height = 10, units = "in", dpi = 300)
dev.off()

# volcano plot
data_for_volanco <- df %>% nest(data = -cell_type) %>% 
    mutate(test_df = map(data, ~.x %>% pivot_transcript())) %>% 
    select(cell_type,test_df) %>% unnest(test_df)

# read all ppcseq results
ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/ppcseq_data/"
ppcseq_files = list.files(path = glue("{ppcseq_dir}"))

ppcseq_df <- tibble(filename = ppcseq_files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(ppcseq_dir, .))) # a new data column
    ) 
ppcseq_df = ppcseq_df %>% unnest(cols = file_contents) %>% select(filename, transcript, tot_deleterious_outliers)
colnames(ppcseq_df)[1] <- "cell_type"

# change cell_type from ....rds to ""
ppcseq_df$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df$cell_type)
# ppcseq_df %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/all_genes_with_outliers.rds")

# combine two dfs
data_for_volanco = data_for_volanco %>% left_join(ppcseq_df) 
data_for_volanco$tot_deleterious_outliers[is.na(data_for_volanco$tot_deleterious_outliers)] <- 0 # change all NA to 0 
data_for_volanco %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/final_sum_table/dataframe_for_volcano.rds")


# filter out the essential transcript (top 5 Pvalue) with outliers
df_label_outliers = data_for_volanco %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    nest(data = -cell_type) %>% 
    mutate(outlier_transcript = map(data, ~ .x %>% arrange(.$PValue))) %>% # 从小到大
    mutate(outlier_transcript = map(outlier_transcript, ~ .x %>% head(5))) %>% 
    select(cell_type, outlier_transcript) %>% 
    unnest(cols = outlier_transcript) %>% 
    select(cell_type, transcript) # select only cell_type and transcript for later left_join with df_for_volcano
# change column name of transcript to symbol (later to label out)



df_label_outliers <- df_label_outliers %>% 
    mutate(symbol = transcript)

# combine two df
data_for_volanco <- data_for_volanco %>% left_join(df_label_outliers, by = c("cell_type", "transcript"))
# replace NA value to ""
data_for_volanco <- data_for_volanco %>% dplyr::mutate(symbol = replace_na(symbol, ""))

dev.off()
dev.new(width = 3000, height = 2000, unit = "px") 
dev.new()
volcano_plot <- data_for_volanco %>%
    ggplot(aes(x = logFC, y = PValue, label = symbol)) +
    geom_point(aes(color = is_significant, size = is_significant, alpha = is_significant)) +
    geom_text_repel(max.overlaps = 30, size = 2) +
    # geom_point(data = subset(df_for_volcano, df_for_volcano$tot_deleterious_outliers > 0),color = "red") +
    # geom_text_repel(aes(logFC, PValue, label = symbol),
    #                 # mapping = aes(label = transcript),
    #                 # color = "red", size = 3,
    #                 box.padding = unit(0.5, "lines"), #字到点的距离
    #                 point.padding = unit(0.8, "lines") #短线段可以省略 +
    #                 segment.color = "red")
    # ) 
    facet_wrap(~ cell_type) +
    custom_theme +
    scale_y_continuous(trans = "log10_reverse") +
    scale_color_manual(values = c("black","#e11f28")) +
    scale_size_discrete(range = c(0, 2))

volcano_plot %>% ggsave(filename = "Volcano_plot.pdf", 
                    path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                    width = 14, height = 10, units = "in", dpi = 300)

volcano_plot %>% ggsave(filename = "Volcano_plot.png", 
                    path = "/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/plot/",
                    width = 14, height = 10, units = "in", dpi = 300)

dev.off()