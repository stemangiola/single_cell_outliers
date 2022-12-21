# Generate 5 summary plots
# aggregate data from all cell types together

library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ggrepel)
library(ggallin)
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
                text = element_text(size = 15),
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

in_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/de_data/" # path to data
out_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/"
    
files = list.files(path = glue("{in_dir}")) # get the file name

df <- files %>% 
    map(~ readRDS(file.path(in_dir, .))) %>% 
    # purrr::reduce(rbind)
    dplyr::bind_rows()

df %>% saveRDS(glue("{out_dir}all_de.rds"))
# df <- readRDS(glue("{out_dir}all_de.rds"))

# generate faced density plot for both 3 methods
# method = c("deseq2", "edgeR_quasi_likelihood","edger_robust_likelihood_ratio")
method = "edgeR_quasi_likelihood"
FDR_column_name = case_when(
    method == "deseq2" ~ "deseq2_padj",
    method == "edgeR_quasi_likelihood" ~ "edgerQLT_FDR",
    method == "edger_robust_likelihood_ratio" ~ "edgerRobust_FDR")
pvalue_column_name = case_when(
    method == "deseq2" ~ "deseq2_pvalue",
    method == "edgeR_quasi_likelihood" ~ "edgerQLT_PValue",
    method == "edger_robust_likelihood_ratio" ~ "edgerRobust_PValue")
logFC_column_name = case_when(
    method == "deseq2" ~ "deseq2_log2FoldChange",
    method == "edgeR_quasi_likelihood" ~ "edgerQLT_logFC",
    method == "edger_robust_likelihood_ratio" ~ "edgerRobust_logFC")

df <- df %>% filter(method == !!method)


df %>% 
    ggplot(aes(x = abundance_RNA_scaled, color = sample)) +
    geom_density() +
    facet_wrap(~ cell_type) +
    scale_x_log10() +
    ggtitle("Density plot of 24 cell types") +
    custom_theme 

# PCA plot
# df1 = df %>% filter(method == "deseq2") %>% 
df1 = df %>% 
    nest(data = -cell_type) %>% 
    mutate(test_df = map(data, ~.x %>% pivot_sample())) %>% 
    select(cell_type,test_df) %>% unnest(test_df)

df1 %>% 
    ggplot(aes(x = PC1, y = PC2, colour = response)) +
    geom_point() +
    facet_wrap(~ cell_type) +
    ggtitle(glue("PCA plot of 24 cell types for method: {method}")) +
    custom_theme

# bar plot
df <- df %>% 
    mutate(abundance_RNA = abundance_RNA %>% as.integer) %>% 
    mutate(is_significant = !!sym(FDR_column_name) < 0.05) #"edgerQLT_FDR"

new_df = df %>% nest(data = -cell_type) %>% 
    mutate(tot_transcript = map(data, ~ .x$transcript %>% n_distinct()))
new_df = new_df %>% select(-data) %>% unnest(tot_transcript)
# sum_table <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/sum_table.rds")
sum_table <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/merged_sum_table.rds")

sum_table <- sum_table %>% filter(method == !!method)

sum_table = new_df %>% left_join(sum_table) 
sum_table <- sum_table %>% mutate(method = !!method)
sum_table[is.na(sum_table)] <- 0
names(sum_table) <- c("cell_type","transcript", "genes_with_outliers","de_genes","outliers_in_top_PValue", "method")
sum_table = sum_table %>% select(cell_type, transcript, de_genes, genes_with_outliers, outliers_in_top_PValue, method)
sum_table <- sum_table %>% pivot_longer(cols = c(transcript, de_genes, genes_with_outliers, outliers_in_top_PValue),
                                        names_to = "variable",
                                        values_to = "Count")

sum_table %>% saveRDS(glue("{out_dir}{method}_final_sum_table.rds"))
# draw bar plot

sum_table$variable = factor(sum_table$variable, levels = c("transcript","de_genes","genes_with_outliers","outliers_in_top_PValue"))

bar_plot = sum_table %>% ggplot(aes(x = variable, y = Count,fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Count)) +
    scale_y_continuous() +
    custom_theme +
    facet_wrap(~ cell_type) +
    xlab("Variable") +
    ylab("Count") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    # ggtitle("summary table of COVID 19 for method: deseq2")
    ggtitle(glue("Summary table of COVID 19 for method: {method}"))


# volcano plot
df2 = df %>% nest(data = -cell_type) %>% 
    mutate(test_df = map(data, ~.x %>% pivot_transcript())) %>% 
    select(cell_type,test_df) %>% unnest(test_df)

# read all ppcseq results
ppcseq_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/ppcseq_data/"

ppcseq_files = list.files(path = glue("{ppcseq_dir}"))
# cell_names = tools::file_path_sans_ext(ppcseq_files) %>% basename() %>% str_replace_all("_ppcseq","")

ppcseq_df <- tibble(filename = ppcseq_files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(ppcseq_dir, .))) # a new data column
    ) 
ppcseq_df = ppcseq_df %>% unnest(cols = file_contents) %>% select(filename, transcript, tot_deleterious_outliers, method) %>% filter(method == !!method)
colnames(ppcseq_df)[1] <- "cell_type"

# change cell_type from ....rds to ""
ppcseq_df$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df$cell_type)
# ppcseq_df$cell_type <- gsub("deseq2_", "", ppcseq_df$cell_type)
ppcseq_df$cell_type <- gsub(glue("{method}_"), "", ppcseq_df$cell_type)
# ppcseq_df %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/all_genes_with_outliers.rds")

# combine two dfs
df_for_volcano = df2 %>% left_join(ppcseq_df) 
df_for_volcano$tot_deleterious_outliers[is.na(df_for_volcano$tot_deleterious_outliers)] <- 0 # change all NA to 0 
# df_for_volcano %>% saveRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/deseq2_df_with_ppcseq.rds")
df_for_volcano %>% saveRDS(glue("{out_dir}{method}_df_with_ppcseq.rds"))


# filter out the essential transcript (top 5 Pvalue) with outliers
df_label_outliers = df_for_volcano %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    nest(data = -cell_type) %>% 
    # mutate(outlier_transcript = map(data, ~ .x %>% arrange(.$PValue))) %>% # 从小到大
    mutate(outlier_transcript = map(data, ~ .x %>% arrange(!!sym(pvalue_column_name)))) %>% # 从小到大
    mutate(outlier_transcript = map(outlier_transcript, ~ .x %>% head(5))) %>% 
    select(cell_type, outlier_transcript) %>% 
    unnest(cols = outlier_transcript) %>% 
    select(cell_type, transcript) # select only cell_type and transcript for later left_join with df_for_volcano
# change column name of transcript to symbol (later to label out)
df_label_outliers <- df_label_outliers %>% 
    mutate(symbol = transcript)

# combine two df
df_for_volcano <- df_for_volcano %>% left_join(df_label_outliers, by = c("cell_type", "transcript"))
# replace NA value to ""
df_for_volcano <- df_for_volcano %>% dplyr::mutate(symbol = replace_na(symbol, ""))

dev.off()
dev.new(width = 4000, height = 3000, unit = "px") 
df_for_volcano %>%
    # ggplot(aes(x = logFC, y = PValue, label = symbol)) +
    ggplot(aes(x = !!sym(logFC_column_name), y = !!sym(pvalue_column_name), label = symbol)) +
    geom_point(aes(color = is_significant, size = is_significant, alpha = is_significant)) +
    geom_text_repel(max.overlaps = 30, size = 2) +
    facet_wrap(~ cell_type) +
    scale_y_continuous(trans = "log10_reverse") +
    scale_color_manual(values = c("black","#e11f28")) +
    scale_size_discrete(range = c(0, 2)) +
    xlab("Log fold change") +
    ylab("P value") +
    # multipanel_theme +
    custom_theme


## Histogram_fold_change_plot
# ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/ppcseq_data/"
# ppcseq_files = list.files(path = glue("{ppcseq_dir}"))
# cell_names = tools::file_path_sans_ext(ppcseq_files) %>% basename() %>% str_replace_all("_ppcseq","")

ppcseq_df_FC <- tibble(filename = ppcseq_files) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename,          # read files into
                               ~ readRDS(file.path(ppcseq_dir, .))) # a new data column
    ) 

colnames(ppcseq_df_FC)[1] <- "cell_type"
ppcseq_df_FC = ppcseq_df_FC %>% unnest(cols = file_contents) %>% filter(method == !!method) %>% 
    select(cell_type, transcript, sample_wise_data, method) %>% 
    unnest(cols = sample_wise_data) 

# change cell_type from ....rds to ""
ppcseq_df_FC$cell_type <- gsub("_ppcseq.rds", "", ppcseq_df_FC$cell_type)
# ppcseq_df_FC$cell_type <- gsub("deseq2_", "", ppcseq_df_FC$cell_type)
ppcseq_df_FC$cell_type <- gsub(glue("{method}_"), "", ppcseq_df_FC$cell_type)

# generate the df for calculating the FC before and after outlier exclusion
ppcseq_df_FC = ppcseq_df_FC %>% 
    distinct(cell_type, transcript, slope_after_outlier_filtering, slope_before_outlier_filtering, method) %>% 
    mutate(FC = slope_after_outlier_filtering/slope_before_outlier_filtering) 
ppcseq_df_FC %>% 
    # saveRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ma.m/single_cell_outliers/single_cell_database/COVID_19/data/final_sum_table/deseq2_FC_df.rds")
    saveRDS(glue("{out_dir}{method}_FC_df.rds"))

ppcseq_df_FC %>% ggplot(aes(x = FC)) +
    geom_histogram(bins = 30) +
    xlim(-10, 10) +
    scale_x_continuous(trans = pseudolog10_trans, guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ cell_type, scales = "free") + 
    # ggtitle("Histogram of fold changes of differential transcripts with and without outliers for method: deseq2") +
    ggtitle(glue("Histogram of fold changes of differential transcripts with and without outliers for method: {method}")) +
    xlab("ratio of fold change after outlier elimination") + 
    custom_theme
