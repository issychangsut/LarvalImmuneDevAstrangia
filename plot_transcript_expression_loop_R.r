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
    dplyr::select(transcript, spID)
  #annotations_subset <- annotations %>%
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
        geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed") +  # Smooth trend line
        labs(title = paste("Expression of", gene_name), x = "Developmental Stage", y = "Expression Level") +
        theme_classic() +  # White background
        theme(legend.position = "none")  # Hide legend if too many samples
      
      # Save plot using the spID name
      ggsave(filename = file.path("~/Desktop/LarvalDevRI22/graphsformerge", paste0(gene_name, ".png")), 
             plot = p, width = 8, height = 6)
      
      return(p)  # Return the plot in case you want to display it in R
    })
  

