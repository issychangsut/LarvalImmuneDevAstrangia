##okay the following is designed to take your tranacripts from the matching file and graph them over time. 
#It includes a loop so that it will go through every trinityID in one file and plot them, 1 for each transcript
#I also wanted to label them with the spID instead of trinityID so that I didnt have to look it up later 

plot_transcript_expression <- function(expression_data, metadata, timepoint_col, timepoint_order, output_dir = "plots") 
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Create output directory if it doesn't exist
  if (!dir.exists("~/Desktop/LarvalDevRI22/graphsformerge")) {
    dir.create("~/Desktop/LarvalDevRI22/graphsformerge")
  }
  
  #load metadata
  metadata <- read.csv("~/Desktop/LarvalDevRI22/LarvalImmunityMeta.csv")
  expression_data<- read.csv("~/Desktop/LarvalDevRI22/matching_transcriptsreads.csv")
  
  # Gather expression data into long format
  long_expr <- expression_data %>%
    pivot_longer(cols = X38:X84, names_to = "SequenceID", values_to = "Expression") 
  
  #my columns in long_expr dont match the file I want to merge with - one is called 
  #sampleID the other is sequenceID so going to manually change them to match in excel
  long_expr<- read.csv("~/Desktop/LarvalDevRI22/long_expr.csv")
  
  # Merge with metadata
  merged_data <- long_expr %>%
    left_join(metadata, by = "SequenceID")
  #I also want the spIDs so that I can name my plots easier 
  annotations_subset <- annotations %>%
    select(transcript, spID)  # Keep only transcript and spID
  
  merged_data <- merged_data %>%
    left_join(annotations_subset, by = "transcript") 
  
  # Manually order timepoints
  merged_data[["developmental_stage"]] <- factor(merged_data[["developmental_stage"]], 
                                               levels = c("zygote", "2_cell", "4_cell", "blastula", 
                                                          "gastrula", "planula_48h_UF", 
                                                          "planula_72h_UF", "planula_96h_UF"))
  
  
  #kept having extra stuff that said NA since i didnt select it so i got rid of all the NAs
  merged_data <- merged_data %>%
    filter(!is.na(developmental_stage))

  # Generate and save plots for each transcript
  merged_data %>%
    split(.$transcript) %>%
    lapply(function(df) {
      gene_name <- ifelse(!is.na(df$spID[1]), df$spID[1], df$transcript[1])  # Use spID if available
      
      p <- ggplot(df, aes(x = developmental_stage, y = Expression, group = SequenceID)) +
        geom_line(aes(color = SequenceID), alpha = 0.6) +  # Line plot
        geom_point(aes(color = SequenceID)) +  # Scatter points
        geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", linetype = "dashed") +  # Smooth trend line
        labs(title = paste("Expression of", gene_name), x = "Developmental Stage", y = "Expression Level") +
        theme_classic() +  # White background
        theme(legend.position = "none")  # Hide legend if too many samples
      
      # Save plot using the spID name
      ggsave(filename = file.path("~/Desktop/LarvalDevRI22/graphsformerge", paste0(gene_name, ".png")), 
             plot = p, width = 8, height = 6)
      
      return(p)  # Return the plot in case you want to display it in R
    })
  
##playing around here
#using a Dynamic Time Warping analysis
  #okay but now i kinda wanna see if any of the genes cluster 
  # Load required libraries
  #install.packages("dtw")
  library(dtw)
  library(pheatmap)
  library(tidyverse)
  
  # Load expression data
  expression_data <- read.csv("~/Desktop/LarvalDevRI22/matching_transcriptsreads.csv", row.names = 1)
  
  # Select only the expression columns (assuming they are numeric and ordered by timepoint)
  expr_matrix <- expression_data %>%
    dplyr::select(X38:X84) %>%
    as.matrix()
  
  # Compute DTW distance matrix
  dtw_dist <- dist(expr_matrix, method = "DTW")
  
  # Perform hierarchical clustering
  hc <- hclust(dtw_dist, method = "complete")
  
  # Plot dendrogram
  plot(hc, main = "Hierarchical Clustering of Transcripts (DTW)")
  
  # Assign clusters (adjust k as needed)
  k <- 5
  clusters <- cutree(hc, k)
  
  # Add cluster info to the data
  cluster_assignments <- data.frame(transcript = rownames(expr_matrix), cluster = clusters)
  
  # Save cluster assignments
  write.csv(cluster_assignments, "~/Desktop/LarvalDevRI22/transcript_clusters_DTW_415k5.csv", row.names = FALSE)
  
  # Heatmap visualization of DTW distances
  #pheatmap(as.matrix(dtw_dist), clustering_method = "complete", main = "DTW Distance Heatmap")
  
##Visulaizing each cluster 
  library(tidyverse)  # Ensure tidyverse is loaded
  
  # Load expression data and cluster assignments
  expression_data <- read.csv("~/Desktop/LarvalDevRI22/matching_transcriptsreads.csv", row.names = 1)
  cluster_assignments <- read.csv("~/Desktop/LarvalDevRI22/transcript_clusters_DTW_415k5.csv")
  
  # Ensure transcript column exists
  expression_data <- expression_data %>%
    rownames_to_column(var = "transcript") %>%
    inner_join(cluster_assignments, by = "transcript")
  
  # Pivot the data from wide to long format
  long_expr <- expression_data %>%
    pivot_longer(cols = starts_with("X"), names_to = "SequenceID", values_to = "expression")
  
  # Remove the "X" from SequenceID in long_expr to match metadata
  long_expr <- long_expr %>%
    mutate(SequenceID = gsub("^X", "", SequenceID))  # Remove leading "X"
  
  # Check unique SequenceID values before merging
  print("Unique SequenceIDs in long_expr:")
  print(unique(long_expr$SequenceID))
  
  # Load metadata
  metadata <- read.csv("~/Desktop/LarvalDevRI22/LarvalImmunityMeta.csv")
  
  # Ensure SequenceID is the same type in both datasets
  metadata$SequenceID <- as.character(metadata$SequenceID)
  long_expr$SequenceID <- as.character(long_expr$SequenceID)
  
  # Join with metadata to get developmental_stage
  long_expr <- long_expr %>%
    left_join(metadata, by = "SequenceID")
  # Manually reorder the developmental stages
  long_expr$developmental_stage <- factor(long_expr$developmental_stage, 
                                          levels = c("zygote", "2_cell", "4_cell", "blastula", 
                                                     "gastrula", "planula_48h_UF", "planula_72h_UF", 
                                                     "planula_96h_UF"))
  
  # Calculate mean expression per cluster and developmental stage
  mean_expr <- long_expr %>%
    group_by(cluster, developmental_stage) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop")
  
  
  # Ensure factor ordering
  mean_expr$developmental_stage <- factor(mean_expr$developmental_stage, 
                                          levels = c("zygote", "2_cell", "4_cell", "blastula", "gastrula", "planula_48h_UF", "planula_72h_UF", 
                                                     "planula_96h_UF"))
  
  # Remove potential NA values
  mean_expr <- mean_expr %>% filter(!is.na(developmental_stage))
  
  # Plot the mean expression for each cluster with a best fit line (LOESS) and confidence interval
  ggplot(mean_expr, aes(x = developmental_stage, y = mean_expression, color = factor(cluster))) +
    geom_line(aes(group = cluster), alpha = 0.8) +  # Line plot for each cluster
    geom_point(aes(shape = factor(cluster)), size = 3) +  # Points for each cluster
    geom_smooth(method = "loess", aes(group = 1), se = TRUE, color = "black", linetype = "dashed") +  # Trendline per facet
    facet_wrap(~ cluster, scales = "free_y") +  # Facet by cluster
    labs(title = "Mean Gene Expression Across Developmental Stages",
         x = "Developmental Stage", y = "Mean Expression Level", color = "Cluster", shape = "Cluster") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
  
   
  
  ##also messing around and trying to look at groupings-- i hate heatmaps these sucked 
  install.packages("pheatmap")
  library(pheatmap)
  
  # Reorder rows based on cluster assignment
  expr_matrix_ordered <- expr_matrix[order(cluster_assignments$cluster), ]
  
  # Create heatmap
  pheatmap(expr_matrix_ordered, scale = "row", clustering_method = "ward.D2")
  
  ##alright now what i want to do is filter to only significant groupings then run gomwu to see which go annotations define each cluster
  mean_expr <- mean_expr %>%
    group_by(cluster) %>%
    mutate(z_score = (mean_expression - mean(mean_expression)) / sd(mean_expression))
  
  z_threshold <- quantile(abs(mean_expr$z_score), 0.90)  # Get top 10%
  significant_clusters <- mean_expr %>%
    filter(abs(z_score) >= z_threshold)

  
  ##Do i need this? need help w this part
  #####i wanna write a for loop saying for each cluster,  count the number of times each go term appears
  library(dplyr)
  GO_enrichment <- 
  
    
    
    go_term_counts <- significant_clusters %>%
    group_by(cluster, GO_term) %>%
    summarise(count = n()) %>%
    arrange(cluster, desc(count))  # Sort by most frequent to see which overall process is mediated by mzt?
 
  
  
  