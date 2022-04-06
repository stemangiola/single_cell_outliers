B_cells_filtered <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/filtered_data/B_cells_filtered.rds")
B_cells_filtered %>% group_by(sample) %>% 
    arrange(median(abundance_RNA))
    
B_cells_filtered %>% with_groups(sample, ~ mutate(.x, median(abundance_RNA))) %>% 
    select(sample, transcript, abundance_RNA, `median(abundance_RNA)`)


B_cells_filtered %>% with_groups(sample, ~ mutate(.x, min(abundance_RNA))) %>% 
    select(sample, transcript, abundance_RNA, `min(abundance_RNA)`)

B_test <- B_cells_filtered %>% keep_abundant(
    .sample = sample,
    .transcript = transcript,
    .abundance = abundance_RNA,
    factor_of_interest = condition,
    minimum_counts = 125,
    minimum_proportion =  0.7) %>% # factor of interest = condition (response/non-response)
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance_RNA) 

B_test %>% with_groups(sample, ~ mutate(.x, median(abundance_RNA))) %>% 
    select(sample, transcript, abundance_RNA, `median(abundance_RNA)`)

# 1. Counts the distinct transcripts
B_test %>% distinct(transcript) # 33 transcripts


# 2. Draw densities for scaled counts
B_test %>%
    # Reshaping
    pivot_longer(cols = c("abundance_RNA", "abundance_RNA_scaled"), names_to = "source", values_to = "abundance") %>%
    
    # Plotting
    ggplot(aes(x = abundance + 1, color = sample)) +
    geom_density() +
    facet_wrap(~source) +
    scale_x_log10() +
    custom_theme

# 3. Histogram of multipliers

B_test %>% ggplot(aes(x = multiplier)) +
    geom_histogram(bins = 30) +
    ggtitle("Histogram of Multiplier for samples")

# 4. Draw PCA

# PCA plot
# Get principal components
B_test_PCA <-
    B_test %>%
    reduce_dimensions(method = "PCA")

B_test_PCA %>%
    pivot_sample() %>%
    ggplot(aes(x = PC1, y = PC2, colour = log(size))) +
    scale_fill_distiller(palette = "Spectral") +
    geom_point() 





PCA_plot <- df_for_PCA %>% 
    ggplot(aes(x = PC1, y = PC2, colour = log(size))) +
    scale_fill_distiller(palette = "Spectral") +
    geom_point() +
    # geom_text_repel(aes(label = sample), show.legend = FALSE) +
    facet_wrap(~ cell_type) +
    ggtitle("PCA plot for 11 cell types for dataset (GSE120575)") 

