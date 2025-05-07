
setwd("~/Desktop/LarvalDevRI22")
# Load necessary libraries
library(tibble)
library(dplyr)
library(tidyr)
library(topGO)
library(ggplot2)

# Read the clustered transcript file (containing transcripts and their assigned cluster)
clusters <- read.csv("~/Desktop/LarvalDevRI22/transcript_clusters_DTW_alltranscripts.csv", header = TRUE, stringsAsFactors = FALSE)
# Read the annotations file (containing transcripts and GO terms)
annos <- read.csv("~/Desktop/LarvalDevRI22/astrangia_annos_clean.csv", header = TRUE, stringsAsFactors = FALSE)

# Clean and prepare data
annos <- subset(annos, select = c("transcript", "GO"))
merged_data <- left_join(clusters, annos, by = "transcript")

merged_data_clean <- merged_data %>%
  separate_rows(GO, sep = ";") %>%
  mutate(GO = trimws(GO))

# Save cleaned merged data
write.csv(merged_data_clean, "~/Desktop/LarvalDevRI22/clustered_transcripts_with_GO_alldata.csv", row.names = FALSE)

# Create gene-to-GO list
geneID2GO <- merged_data_clean %>%
  group_by(transcript) %>%
  summarise(GO = list(unique(GO))) %>%
  deframe()

# Make sure GO terms are cleaned properly
geneID2GO_cleaned <- lapply(geneID2GO, function(go_terms) {
  go_terms <- unlist(strsplit(go_terms, ";"))
  go_terms <- trimws(go_terms)
  return(go_terms)
})

# Define all transcripts
all_genes <- unique(merged_data_clean$transcript)

# Function to run enrichment for a specific cluster
run_GO_enrichment <- function(cluster_id) {
  cat("Processing Cluster", cluster_id, "...\n")
  
  # Get genes in the selected cluster
  interesting_genes <- merged_data_clean %>%
    filter(cluster == cluster_id) %>%
    pull(transcript) %>%
    unique()
  
  # Create named gene list
  geneList <- factor(as.integer(all_genes %in% interesting_genes))
  names(geneList) <- all_genes
  
  # Check gene-GO overlap
  overlap <- length(intersect(names(geneList), names(geneID2GO_cleaned)))
  cat("  Overlapping genes with GO terms:", overlap, "\n")
  if (overlap == 0) {
    cat("  Skipping cluster due to no overlap.\n")
    return(NULL)
  }
  
  # Create topGO object
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x) x == 1,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO_cleaned)
  
  # Run enrichment test
  result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Extract results
  enrichment_results <- GenTable(GOdata, classicFisher = result, topNodes = 500)
  
  # Save results
  output_file <- paste0("~/Desktop/LarvalDevRI22/GO_enrichment_all_data_Cluster_", cluster_id, ".csv")
  write.csv(enrichment_results, output_file, row.names = FALSE, quote = TRUE)
  
  # Plotting
  plot_GO_terms(cluster_id)
}

# Plotting function
plot_GO_terms <- function(cluster) {
  file_path <- paste0("~/Desktop/LarvalDevRI22/GO_enrichment_all_data_Cluster_", cluster, ".csv")
  
  if (file.exists(file_path)) {
    go_data <- read.csv(file_path)
    go_data <- go_data[go_data$classicFisher != "< 1e-30", ]  # Optional: remove extreme p-values that can't be converted
    
    if (nrow(go_data) == 0) {
      message("No significant GO terms for Cluster: ", cluster)
      return(NULL)
    }
    
    go_data$classicFisher <- as.numeric(go_data$classicFisher)
    
    p <- ggplot(go_data, aes(x = reorder(Term, -log10(classicFisher)), y = -log10(classicFisher))) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = paste("Top GO Terms for Cluster", cluster),
           x = "GO Term",
           y = "-log10(p-value)") +
      theme_minimal()
    
    ggsave(filename = paste0("~/Desktop/LarvalDevRI22/GO_Plot_all_data_Cluster_", cluster, ".png"), plot = p, width = 8, height = 6)
    
    return(p)
  } else {
    message("File not found for Cluster: ", cluster)
    return(NULL)
  }
}

# Run for all unique clusters
all_clusters <- sort(unique(merged_data_clean$cluster))
for (cid in all_clusters) {
  try(run_GO_enrichment(cid))
}
###########
###########Alrigthy now im looking at the specific immune ones 


# Load necessary libraries
library(dplyr)
library(ggplot2)

# Directory with your enrichment files
output_dir <- "/Users/isabellachangsut/Desktop/LarvalDevRI22"

# Define immune-related keywords for flexible matching (including partial words)
immune_keywords <- c("immun", "defens", "cytokin", "antigen", "inflam", "infect", "innate", "autophagy", "ubiquitin", "apopt", "autoph", "antibact", "NFKB","patho", "virus", "viral", "interferon", "lipopol", "LPS", "stress")

# Get all enrichment CSVs in the output folder
enrichment_files <- list.files(path = output_dir, pattern = "^GO_enrichment_all_data_Cluster_\\d+\\.csv$", full.names = TRUE)

# Function to filter and plot immune-related GO terms
plot_immune_GO_terms <- function(file_pattern) {
  # Extract cluster ID from file name
  cluster_id <- gsub(".*Cluster_(\\d+)\\.csv$", "\\1", file_path)
  
  # Read the enrichment data
  go_data <- read.csv("~/Desktop/LarvalDevRI22/GO_enrichment_Cluster(\\d+)\\.csv$")
  
  # Clean p-values (remove "< 1e-30" or similar)
  go_data <- go_data[go_data$classicFisher != "< 1e-30", ]
  go_data$classicFisher <- as.numeric(go_data$classicFisher)
  
# Filter for immune-related GO terms using the more inclusive keywords, if there are no terms then say no immune related go terms for the cluster #
  immune_terms <- go_data %>%
    filter(grepl(paste(immune_keywords, collapse = "|"), Term, ignore.case = TRUE))
  
  if (nrow(immune_terms) == 0) {
    message("No immune-related GO terms for Cluster: ", cluster_id)
    return(NULL)
  }
  
  # Create plot for immune-related terms
  p <- ggplot(immune_terms, aes(x = reorder(Term, -log10(classicFisher)), y = -log10(classicFisher))) +
    geom_bar(stat = "identity", fill = "black") +
    coord_flip() +
    labs(title = paste("Immune-Related GO Terms - Cluster", cluster_id),
         x = "GO Term", y = "-log10(p-value)") +
    theme_minimal()
  
  # Save the plot as PNG
  plot_file <- file.path(output_dir, paste0("Immune_GO_Plot_Cluster_", cluster_id, ".png"))
  ggsave(plot_file, plot = p, width = 8, height = 6)
  
  return(p)
}

# Loop through each enrichment file and generate immune-related plots
lapply(enrichment_files, plot_immune_GO_terms)
pdf(file.path(output_dir, "All_Immune_GO_Plots.pdf"), width = 8, height = 6)







#################################################################
setwd("~/Desktop/LarvalDevRI22")
# Load necessary libraries
library(dplyr)
library(tidyr)
library(topGO)
# Read the clustered transcript file (containing transcripts and their assigned cluster)
clusters <- read.csv("~/Desktop/LarvalDevRI22/transcript_clusters_DTW_alltranscripts.csv", header = TRUE, stringsAsFactors = FALSE)

# Read the annotations file (containing transcripts and GO terms)
annos <- read.csv("~/Desktop/LarvalDevRI22/astrangia_annos_clean.csv", header = TRUE, stringsAsFactors = FALSE)

# Keep only relevant columns: transcript and GO term
annos <- annos %>%
  dplyr::select(transcript, GO)

# Merge clusters with GO annotations
merged_data <- left_join(clusters, annos, by = "transcript")

# Save the merged file for reference
write.csv(merged_data, "clustered_transcripts_with_GO.csv", row.names = FALSE)


# Load necessary libraries
library(tibble)
library(dplyr)
library(tidyr)
library(topGO)

# Subset to keep only necessary columns
annos <- subset(annos, select = c("transcript", "GO"))

# Merge annotation with cluster assignments
merged_data <- left_join(clusters, annos, by = "transcript")

# Split multiple GO terms into separate rows (separated by ;)
merged_data_clean <- merged_data %>%
  separate_rows(GO, sep = ";") %>%
  mutate(GO = trimws(GO))  # Trim any whitespace

# Save for reference
write.csv(merged_data_clean, "~/Desktop/LarvalDevRI22/clustered_transcripts_with_GO_alldata.csv", row.names = FALSE)

# Count how many transcripts have/don't have GO terms
cat("Transcripts with GO terms:", sum(!is.na(merged_data$GO)), "\n")
cat("Transcripts without GO terms:", sum(is.na(merged_data$GO)), "\n")

# Count the number of times each GO term appears per cluster
go_counts <- merged_data %>%
  group_by(cluster, GO) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(cluster, desc(count))

write.csv(go_counts, "GO_terms_per_cluster_alldatak2.csv", row.names = FALSE)

# Create a gene-to-GO list
geneID2GO <- merged_data_clean %>%
  group_by(transcript) %>%
  summarise(GO = list(unique(GO))) %>%
  deframe()

# OPTIONAL (extra clean-up if some GO terms are still joined by ";")
geneID2GO_cleaned <- lapply(geneID2GO, function(go_terms) {
  go_terms <- unlist(strsplit(go_terms, ";"))
  go_terms <- trimws(go_terms)
  return(go_terms)
})

# Pick a cluster to analyze for enrichment (e.g., cluster 1)
selected_cluster <- 2

# Get genes in the selected cluster
interesting_genes <- merged_data_clean %>%
  filter(cluster == selected_cluster) %>%
  pull(transcript) %>%
  unique()

# Define all genes (background)
all_genes <- unique(merged_data_clean$transcript)

# Create named binary vector (1 = in cluster, 0 = not)
geneList <- factor(as.integer(all_genes %in% interesting_genes))
names(geneList) <- all_genes

# Check overlap
cat("Overlap between geneList and geneID2GO:", length(intersect(names(geneList), names(geneID2GO_cleaned))), "\n")

# Create the topGO object
GOdata <- new("topGOdata",
              ontology = "BP",  # Change to "MF" or "CC" if desired
              allGenes = geneList,
              geneSelectionFun = function(x) x == 1,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO_cleaned)

# View summary
GOdata

# Run GO enrichment analysis
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Summarize the results
summary(result)

# Save results
output_file <- paste0("GO_enrichment_Cluster_", cluster_id, ".csv")
write.csv(enrichment_results, output_file, row.names = FALSE, quote = TRUE)


# Plotting
plot_GO_terms(cluster_id)
}

# Plotting function
plot_GO_terms <- function(cluster) {
  file_path <- paste0("GO_enrichment_Cluster_", cluster, ".csv")
  
  if (file.exists(file_path)) {
    go_data <- read.csv(file_path)
    go_data <- go_data[go_data$classicFisher != "< 1e-30", ]  # Optional: remove extreme p-values that can't be converted
    
    if (nrow(go_data) == 0) {
      message("No significant GO terms for Cluster: ", cluster)
      return(NULL)
    }
    
    go_data$classicFisher <- as.numeric(go_data$classicFisher)
    
    p <- ggplot(go_data, aes(x = reorder(Term, -log10(classicFisher)), y = -log10(classicFisher))) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = paste("Top GO Terms for Cluster", cluster),
           x = "GO Term",
           y = "-log10(p-value)") +
      theme_minimal()
    
    ggsave(filename = paste0("GO_Plot_Cluster_", cluster, ".png"), plot = p, width = 8, height = 6)
    
    return(p)
  } else {
    message("File not found for Cluster: ", cluster)
    return(NULL)
  }
}

# Run for all unique clusters
all_clusters <- sort(unique(merged_data_clean$cluster))
for (cid in all_clusters) {
  try(run_GO_enrichment(cid))
}



