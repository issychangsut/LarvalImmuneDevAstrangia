# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Set your mf working directory but also i used complete paths anyway bc my files are a mess <3
setwd("~/Desktop/LarvalDevRI22")

# Load your cluster assignments
clusters <- read.csv("~/Desktop/LarvalDevRI22/transcript_clusters_DTW_alltranscripts.csv", header = TRUE, stringsAsFactors = FALSE)
# Load your annotations with GO terms
annos <- read.csv("~/Desktop/LarvalDevRI22/astrangia_annos_clean.csv", header = TRUE, stringsAsFactors = FALSE)
annos <- subset(annos, select = c("transcript", "anno"))
annos <- as_tibble(annos)
# Load your original expression matrix used for clustering (your normalized filtered reads)
expression_matrix <- read.csv("~/Desktop/LarvalDevRI22/larvalnormalizedfilteredreads.csv", header = TRUE)
# Join clusters with GO annotations
merged_data <- left_join(clusters, annos, by = "transcript")

# Define immune-related keywords so that you can break the clusters into immune and non immune
immune_keywords <- c("immun", "defens", "cytokin", "antigen", "inflam", "infect", "innate", 
                     "autophagy", "ubiquitin", "apopt", "autoph", "antibact", "NFKB",
                     "patho", "virus", "viral", "interferon", "lipopol", "LPS", "stress")

#labeling any of them as immune transcripts
immune_transcripts <- annos %>%
  filter(grepl(paste(immune_keywords, collapse = "|"), anno, ignore.case = TRUE)) %>%
  pull(transcript) %>%   # This matches the 'transcript' field in your other files
  unique()

#make a matrix with expression and the cluster assignment
expr_clustered <- expression_matrix %>%
  inner_join(clusters, by = "transcript")
# only need the transcript and cluster assignment
feature_cols <- setdiff(names(expr_clustered), c("transcript", "cluster"))

####
####
### Computing the centroids for each cluster
####
####

centroids <- expr_clustered %>%
  group_by(cluster) %>%
  summarise(across(all_of(feature_cols), mean), .groups = "drop")

####
####
### Computing the Euclidean distance from centroids to each transcript (row)
### this will let me know which are the top drivers of the patterns we see (lower distance = stronger driver)
####
####

expr_clustered_with_dist <- expr_clustered %>%
  rowwise() %>%
  mutate(distance = sqrt(sum((c_across(all_of(feature_cols)) - 
                                centroids[centroids$cluster == cluster, feature_cols])^2))) %>%
  ungroup()

# Adding the group labels (immune vs other)
expr_clustered_with_dist <- expr_clustered_with_dist %>%
  mutate(type = ifelse(transcript %in% immune_transcripts, "Immune", "Other"))


####
####
### Plotting onto a violin plot, one per cluster
####
####

# Making sure cluster is a factor so ggplot facets it correctly
violin_data <- expr_clustered_with_dist %>%
  filter(cluster %in% c(1, 2)) %>%
  mutate(cluster = factor(cluster),
         type = factor(type, levels = c("Immune", "Other")))

cluster_1_data <- expr_clustered_with_dist %>% filter(cluster == 1)
# Create violin plot for Cluster 1 (showing immune vs other)
p_cluster_1 <- ggplot(cluster_1_data, aes(x = type, y = distance, fill = type)) +
  geom_violin(trim = FALSE, scale = "width", color = "black") + stat_summary(fun = "mean", geom = "crossbar", color = "black") +
  labs(title = "Distance from Cluster Centroid: Immune vs Other Genes (Cluster 1)",
       x = "Gene Type", y = "Distance from Centroid") +
  scale_fill_manual(values = c("Immune" = "#FFCCCB", "Other" = "aliceblue")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Print the plot
print(p_cluster_1)

# save the plot if it looks good
ggsave("~/Desktop/Cluster1_Immune_vs_Other_Violin_Plot.png", plot = p_cluster_1, width = 5, height = 8)

#####
#####
### SECOND CLUSTER NOW
#####
#####

# Subset data for cluster 2
cluster_2_data <- expr_clustered_with_dist %>% filter(cluster == 2)

# Add immune label based on matching keywords in 'anno' (annotation column)
cluster_2_data <- cluster_2_data %>%
  mutate(type = ifelse(grepl(paste(immune_keywords, collapse = "|"), anno, ignore.case = TRUE), 
                       "Immune", "Other"))

# Plot the violin plot for cluster 2
p_cluster_2 <- ggplot(cluster_2_data, aes(x = type, y = distance, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  labs(title = "Cluster 2: Immune vs. Other Genes",
       x = "Gene Type", y = "Distance from Centroid") + stat_summary(fun = "mean", geom = "crossbar", color = "black") +
  scale_fill_manual(values = c("Immune" = "#FFCCCB", "Other" = "aliceblue")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Print the plot
print(p_cluster_2)

# Save the plot to a file
ggsave("~/Desktop/Cluster2_Immune_vs_Other_Violin.png", plot = p_cluster_2, width = 5, height = 8)


####
####
### FACETING
####
####

# Filter for clusters 1 and 2 and convert cluster to factor
facet_data <- expr_clustered_with_dist %>%
  filter(cluster %in% c(1, 2)) %>%
  mutate(cluster = factor(cluster))

# Create the enhanced violin plot
faceted_bothclusters <- ggplot(facet_data, aes(x = type, y = distance, fill = type)) +
  geom_violin(trim = FALSE, scale = "width", color = "black") + stat_summary(fun = "mean", geom = "crossbar", color = "black") +
  labs(
    title = "Distance from Centroid by Gene Type",
    subtitle = "Immune-related transcripts vs. all other transcripts",
    x = "Gene Type",
    y = "Distance from Centroid"
  ) +
  scale_fill_manual(values = c("Immune" = "#FFCCCB", "Other" = "aliceblue")) +
  scale_color_manual(values = c("Immune" = "#FFCCCB", "Other" = "aliceblue")) +
  facet_wrap(~ cluster, nrow = 1, labeller = label_both) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13)
  )
# Checking if it looks okay
print(faceted_bothclusters)

# Save and print
ggsave("~/Desktop/LarvalDevRI22/BothClusters_ImmuneVsOtherViolin.png", plot = faceted_bothclusters, width = 10, height = 13)


######
####### Count immune-related genes in cluster 2
immune_count_cluster2 <- expr_clustered_with_dist %>%
  filter(cluster == 2, type == "Immune") %>%
  distinct(transcript) %>%
  nrow()

immune_count_cluster2## there are 5

####
####
### Curious about which transcripts are driving the patterns of each cluster
#####
#####

### Cluster 1
# Picking the closest 10% to centroid
threshold_1 <- quantile(cluster_1_data$distance, 0.10)
top_drivers_cluster1 <- cluster_1_data %>%
  filter(distance <= threshold_1)

# Save the top-driving genes
view(top_drivers_cluster1)

### Same for cluster 2 now
threshold_2 <- quantile(cluster_2_data$distance, 0.10)  # 10th percentile
top_drivers_cluster2 <- cluster_2_data %>%
  filter(distance <= threshold_2)
view(top_drivers_cluster2)

# Combine the results for both clusters
top_drivers_combined <- bind_rows(
  top_drivers_cluster1 %>% mutate(cluster = 1),
  top_drivers_cluster2 %>% mutate(cluster = 2))

# If you want to create a violin plot for both clusters, here's an example:
# Create a combined violin plot for both clusters Immune vs Other
violin10 <- ggplot(top_drivers_combined, aes(x = type, y = distance, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(width = 0.05, alpha = 0.7, size = 1) +
  facet_wrap(~cluster) +
  labs(title = "Top Drivers of Clusters 1 and 2: Immune vs Other",
       x = "Gene Type", y = "Distance from Centroid") +
  scale_fill_manual(values = c("Immune" = "#FFCCCB", "Other" = "aliceblue")) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")

# checking to see if it looks okay
print(violin10)

# Save the plot for combined clusters
#ggsave("~/Desktop/LarvalDevRI22/Violin_Top_DriversAllCLusters.png", width = 10, height = 13) its uglieeeee/unimportant ;)

####
####
###okay all done 
####
####