# aggregate data from all cell types together
library(tidyverse)
library(glue)
library(ggplot2)
library(tidybulk)
library(ggrepel)

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
        # panel.grid.major = element_line(size = 0.2),
        # panel.grid.minor = element_line(size = 0.1),
        # text = element_text(size = 12),
        # legend.position = "bottom",
        strip.background = element_blank(),
        # axis.title.x = element_text(margin = margin(
        #   t = 10,
        #   r = 10,
        #   b = 10,
        #   l = 10
        # )),
        # axis.title.y = element_text(margin = margin(
        #   t = 10,
        #   r = 10,
        #   b = 10,
        #   l = 10
        # )),
        # axis.text.x = element_text(
        #   angle = 30,
        #   hjust = 1,
        #   vjust = 1, 
        # )
      )
  )

in_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/de_data/" # path to data
files = list.files(path = glue("{in_dir}")) # get the file name

df <- files %>% 
  map(~ readRDS(file.path(in_dir, .))) %>% 
  reduce(rbind)
# df %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/all_de.rds")

# generate faced density plot
df %>% 
  ggplot(aes(x = abundance_RNA_scaled, color = sample)) +
  geom_density() +
  facet_wrap(~ cell_type) +
  scale_x_log10() +
  ggtitle("Density plot of 24 cell types") +
  custom_theme 

# PCA plot
df1 = df %>% nest(data = -cell_type) %>% 
  mutate(test_df = map(data, ~.x %>% pivot_sample())) %>% 
  select(cell_type,test_df) %>% unnest(test_df)

df1 %>% 
  ggplot(aes(x = PC1, y = PC2, colour = response)) +
  geom_point() +
  facet_wrap(~ cell_type) +
  ggtitle("PCA plot of 24 cell types") +
  custom_theme

# bar plot
df <- df %>% 
  mutate(abundance_RNA = abundance_RNA %>% as.integer) %>% 
  mutate(is_significant = FDR < 0.05) 
  
new_df = df %>% nest(data = -cell_type) %>% 
  mutate(tot_transcript = map(data, ~ .x$transcript %>% n_distinct()))
new_df = new_df %>% select(-data) %>% unnest(tot_transcript)
sum_table <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/sum_table.rds")


sum_table = new_df %>% left_join(sum_table) 
sum_table[is.na(sum_table)] <- 0
names(sum_table) <- c("cell_type","transcript", "genes_with_outliers","de_genes","outliers_in_top_PValue")
sum_table = sum_table %>% select(cell_type, transcript, de_genes, genes_with_outliers, outliers_in_top_PValue)
sum_table <- sum_table %>% pivot_longer(cols = c(transcript, de_genes, genes_with_outliers, outliers_in_top_PValue),
                                        names_to = "variable",
                                        values_to = "Count")
sum_table %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/final_sum_table.rds")
# draw bar plot

sum_table$variable = factor(sum_table$variable, levels = c("transcript","de_genes","genes_with_outliers","outliers_in_top_PValue"))

bar_plot = sum_table %>% ggplot(aes(x = variable, y = Count,fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count)) +
  scale_y_continuous() +
  custom_theme +
  facet_wrap(~ cell_type) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("summary table of COVID 19")


# volcano plot
df2 = df %>% nest(data = -cell_type) %>% 
  mutate(test_df = map(data, ~.x %>% pivot_transcript())) %>% 
  select(cell_type,test_df) %>% unnest(test_df)

# read all ppcseq results
ppcseq_dir = "/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/ppcseq_data/"
ppcseq_files = list.files(path = glue("{ppcseq_dir}"))
cell_names = tools::file_path_sans_ext(ppcseq_files) %>% basename() %>% str_replace_all("_ppcseq","")

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
df_for_volcano = df2 %>% left_join(ppcseq_df) 
df_for_volcano$tot_deleterious_outliers[is.na(df_for_volcano$tot_deleterious_outliers)] <- 0 # change all NA to 0 
df_for_volcano %>% saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/COVID_19/data/final_sum_table/df_with_ppcseq.rds")


# filter out the essential transcript (top 5 Pvalue) with outliers
df_label_outliers = df_for_volcano %>% 
    filter(tot_deleterious_outliers > 0) %>% 
    nest(data = -cell_type) %>% 
    mutate(outlier_transcript = map(data, ~ .x %>% arrange(.$PValue))) %>% # 从小到大
    mutate(outlier_transcript = map(outlier_transcript, ~ .x %>% head(5))) %>% 
    select(cell_type, outlier_transcript) %>% 
    unnest(cols = outlier_transcript)

dev.off()
dev.new(width = 3000, height = 1500, unit = "px") 
df_for_volcano %>%
  ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
  geom_point(data = df_for_volcano) +
  geom_point(data = subset(df_for_volcano, df_for_volcano$tot_deleterious_outliers > 0),color = "red") +
  geom_text_repel(data = df_label_outliers, aes(logFC, PValue, label = transcript),
            # mapping = aes(label = transcript), 
            color = "red", size = 3, 
            box.padding = unit(0.5, "lines"), #字到点的距离
            point.padding = unit(0.8, "lines"), #短线段可以省略
            segment.color = "red"
            # show.legend = F
            ) + # segment.colour = NA, #不显示线段
    facet_wrap(~ cell_type) +
  scale_y_continuous(trans = "log10_reverse") +
    custom_theme


