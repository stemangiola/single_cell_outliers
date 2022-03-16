# investigate the low abundance samples
install.packages("plotly")
library(plotly)
B_cells_tidy <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/tidy_data/B_cells_tidy.rds")
B_cells_tidy <- B_cells_tidy %>% 
    keep_abundant(.sample = sample,
                  .transcript = transcript,
                  .abundance = abundance_RNA,
                  factor_of_interest = condition) %>% # factor of interest = condition (response/non-response)
    scale_abundance(.sample = sample,
                    .transcript = transcript,
                    .abundance = abundance_RNA) 

# find the most abundant genes
# Pick a sample with low transcription. Find the most abundant  gene, tell me what is its count
B_cells_tidy %>% filter(sample == "Pre_P1") %>% arrange(desc(.$abundance_RNA))


# Do for each different cell_types ---> filter out the low abundance samples using makeflow 
# save them into rds and repeat all the process again
df = 
    B_cells_tidy %>% nest(data = -sample) %>% 
    mutate(size = map_dbl(data, ~ .x$abundance_RNA %>% sum())) %>%
    mutate(ratio_from_largest = size/median(size)) %>%
    filter(ratio_from_largest > 0.1) %>% 
    unnest (cols = data) %>%
    pivot_longer(cols = c("abundance_RNA", "abundance_RNA_scaled"), names_to = "source", values_to = "abundance") %>%
    
    # Plotting
    ggplot(aes(x = abundance +1, color = sample)) +
    geom_density() +
    facet_wrap(~source) +
    scale_x_log10() +
    custom_theme
    

B_cells_tidy

df <- df %>% arrange(.$size)

plot2 %>% plotly::ggplotly()

ggplotly(p1)
