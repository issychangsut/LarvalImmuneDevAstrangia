####First we are going to set the wd 
setwd("~/Desktop/LarvalDevRI22")

###Downloading and loading packages 
#install("Biostrings", force = TRUE)
#install("DESeq2", force = TRUE)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(tidyr)
library(DescTools)
library(data.table)

###############################################################################
##PREPPING FILES FOR DESEQ
##based on mishas github:
##https://github.com/z0on/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt
###############################################################################


##reading in the counts file which is the read count matrix
reads=read.delim("~/Desktop/LarvalDevRI22/First_tryBenthics 2024/larvalallcounts.txt", header=TRUE, sep = "\t")
colnames(reads)=gsub(".trim.sam.counts", '', colnames(reads), fixed=TRUE)
colnames(reads) = gsub(pattern = "\\.([[:alpha:]]).*", replacement= "", x=colnames(reads))
reads2=reads[,-1]
rownames(reads2)=reads[,1]

##reading in the metadata from the experiment
#meta=read.csv("LarvalImmunityMeta.csv", header=TRUE, sep=",")
#meta=meta[order(meta$sequenceID),]
#colnames(meta)

##looking at the file and renaming 
#final_reads =reads2 %>% 
 # relocate(order(colnames(.)))
#final_reads=final_reads[,c(1:47)]
#write.csv(meta,"larval_final_meta.csv")

read.csv("~/Desktop/LarvalDevRI22/larval_final_meta.csv")
##Filtering the reads based on mean minimum
means <- apply(final_reads,1,mean)
table(means>3)

means3 <- names(means[means>3])
head(means3)
length(means3)

countFilt <- final_reads[row.names(final_reads) %in% means3,]
head(countFilt)

totalCountsFilt <- colSums(countFilt)
totalCountsFilt

##reporting read count distributions 
min(totalCountsFilt) #0?
max(totalCountsFilt) #619477
mean(totalCountsFilt) #267219

##The following is because I got an error that said every gene had at least one zero so it couldnt do log
##bc of this, i am adding a 1 to each value but want to keep my column and row names
#Here I am convert the matrix to tibble and make the rownames a column so that i can add one to only my reads
count_tibble <- as_tibble(countFilt, rownames = NA) %>% rownames_to_column()%>%
  mutate(across(2:48, ~ . + 1))
##Convert tibble back to matrix then make the first column row names again <3
countFilt2 <- count_tibble%>%column_to_rownames(var = "rowname")%>%as.matrix()

##now noramlizing the reads and getting dispersion estimates 
#dds <- DESeqDataSetFromMatrix(countData = countFilt2, 
#colData = meta,
#design = ~ developmental_stage)

#dds <- estimateSizeFactors(dds)        
#dds <- estimateDispersions(dds)

#vst <- getVarianceStabilizedData(dds)
#write.csv(vst, file = "larvalnormalizedfilteredreads.csv")

##format our uniprot info
annos<-read.delim("astrangia_annos_clean.txt", sep="\t")
annos<-annos[,c(1:3,5)]
colnames(annos)[4]="GO"
colnames(annos)[3]="name"



##modeling for developmental stage = effect of developmental stage 
dds1 <- model.matrix( ~ developmental_stage, meta)
object1 = DESeqDataSetFromMatrix(countData= countFilt2, colData= meta, design = ~ developmental_stage , ignoreRank= TRUE )

object1 <- estimateSizeFactors(object1)

object1 <- estimateDispersions(object1, modelMatrix = dds1)
dds2 <- nbinomWaldTest(object1, maxit=500, modelMatrix = dds1)  
resultsNames(dds2)


#Running contrasts of interest
#
#First, making a function to run the contrasts
StageComp <- function(Stage1, Stage2, data) {
  # Check if Stage1 and Stage2 are character vectors
  if (!is.character(Stage1) || !is.character(Stage2)) {
    stop("Arguments Stage1 and Stage2 must be of type character.")
  }
  
  # Perform differential expression analysis
  DevOneToTwoCell <- results(data, contrast = list(c(Stage1), c(Stage2)))
  DevOneOrdered <- DevOneToTwoCell[order(DevOneToTwoCell$padj),]
  
  # Generate output file name
  output_file <- paste0(Stage1, Stage2, ".csv")
  
  # Write results to CSV file
  write.csv(DevOneOrdered, file = output_file)
}

# Zygote vs Two Cell
DevOneToTwoCell<- results(dds2,contrast=list(c("developmental_stagezygote"), c("Intercept")))
DevOneOrdered = DevOneToTwoCell[order(DevOneToTwoCell$padj),]
write.csv(DevOneOrdered, file = "DESeq_Count_Developmental_stage1.csv")

# Two cell vs 4 cell
StageComp("Intercept","developmental_stage4_cell", dds2)
DevTwo=read.csv("~/Desktop/LarvalDevRI22/Interceptdevelopmental_stage4_cell.csv", header=TRUE)

# Four cell vs blastula
StageComp("developmental_stage4_cell","developmental_stageblastula", dds2)
DevThree=read.csv("developmental_stage4_celldevelopmental_stageblastula.csv", header=TRUE)

#Blastula vs glastula
StageComp("developmental_stageblastula","developmental_stagegastrula", dds2)
DevFour=read.csv("developmental_stageblastuladevelopmental_stagegastrula.csv", header=TRUE)

#Gastrula vs planula 48
StageComp("developmental_stagegastrula","developmental_stageplanula_48h", dds2)
DevFive=read.csv("developmental_stagegastruladevelopmental_stageplanula_48h.csv", header=TRUE)

#Planula 48 vs planula 72
StageComp("developmental_stageplanula_48h", "developmental_stageplanula_72h_UF", dds2)
DevSix=read.csv("developmental_stageplanula_48hdevelopmental_stageplanula_72h_UF.csv", header=TRUE)

#Planula 72 vs 96
StageComp("developmental_stageplanula_72h_UF","developmental_stageplanula_96h_UF", dds2)
DevSeven=read.csv("developmental_stageplanula_72h_UFdevelopmental_stageplanula_96h_UF.csv", header=TRUE)


###############################################################################
########################Looking at specific genes##############################
###############################################################################


##now im doing specific immune shit and reading in my merged data from the code I made "mathcingtranscripts.R" :)
onlyPolr3d<- read.csv("matching_transcriptsreads.csv", row.names=NULL)
# Step 1: Remove the 'X' from the column names
colnames(onlyPolr3d) <- gsub("^X", "", colnames(onlyPolr3d))
# Step 2: If needed, rename columns, especially if the sequence ID is not already in a column
# Example: Rename the first column to sequenceID if it's not already named correctly
colnames(onlyPolr3d)[1] <- "sequenceID"
# Step 1: Correctly filter by sequenceID
subset_data <- subset(onlyPolr3d, sequenceID == "TRINITY_DN10125_c0_g1_i13")

# Step 3: Save the subsetted data to a new CSV
write.csv(subset_data, "subsetonlyPolr3d_.csv", row.names = FALSE)

# Step 4: Check the result
head(subset_data)


# Set the first column as row names and remove it
rownames(subset_data) <- subset_data[, 1]
subset_data <- subset_data[, -1]
normaltransform<- as.data.frame(t(subset_data))
normaltransform
normaltransform$sequenceID<- 38:(37 + nrow(normaltransform))
mergedPolr3d<- merge(normaltransform, meta, by="sequenceID")



###############################################################################
###############################################################################
source("makinggraphsforimmunegene.R")

# Print merged data to check the result
head(mergedPolr3d)



names(mergedPolr3d)[2]<-"Polr3d"
mergedPolr3d<- subset(mergedPolr3d, developmental_stage%in%c("zygote","2_cell","4_cell","blastula","gastrula", "planula_48h", "planula_72h_UF", "planula_96h_UF"))
mergedPolr3d$developmental.stage<- factor(mergedPolr3d$developmental_stage, levels=c("zygote","2_cell","4_cell","blastula","gastrula", "planula_48h", "planula_72h_UF", "planula_96h_UF"),ordered = TRUE)
write.csv(mergedPolr3d, "onlyPolr3dreads.csv")

## Use ggplot to create the plot
Polr3d_plot <- ggplot(data = mergedPolr3d, aes(x = developmental_stage, y = Polr3d)) + 
  geom_smooth(aes(group = 1)) + 
  geom_point() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

# Print the plot
print(Polr3d_plot)
ggsave("~/Desktop/LarvalDevRI22/Figures/Polr3dacrossstagesFED.pdf", Polr3d_plot, height = 4, width = 8)

# Ensure developmental_stage is an ordered factor if it's not already
mergedbby$developmental_stage <- factor(mergedbby$developmental_stage, 
                                        levels = c("zygote", "2_cell", "4_cell", "blastula", "gastrula", 
                                                   "planula_48h", "planula_72h_UF", "planula_96h_UF"))
# Graphing immune genes over time 
#DMBT1
normal<- read.csv("larvalnormalizedfilteredreads.csv", row.names=1)
normalfil<- subset(normal,rownames(normal)%in%c("TRINITY_DN10016_c1_g1_i1"))##row53DMBT1
normaltransform<- as.data.frame(t(normalfil))
normaltransform
normaltransform$sequenceID<- 38:84
mergedbby<- merge(normaltransform, meta, by="sequenceID")

names(mergedbby)[2]<-"DMBT1"
mergedbby<- subset(mergedbby, developmental_stage%in%c("zygote","2_cell","4_cell","blastula","gastrula", "planula_48h", "planula_72h_UF", "planula_96h_UF"))
mergedbby$developmental.stage<- factor(mergedbby$developmental_stage, levels=c("zygote","2_cell","4_cell","blastula","gastrula", "planula_48h", "planula_72h_UF", "planula_96h_UF"),ordered = TRUE)
write.csv(mergedbby, "mergedbby.csv")
marged<-read.csv("mergedbby.csv")
# Create the plot
DMBT1 <- ggplot(data = mergedbby, aes(x = developmental_stage, y = DMBT1)) + 
  geom_smooth(aes(group = 1), method = "loess", se = FALSE) + 
  geom_point() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_blank(),    # Remove external ticks on the y-axis
        axis.text.y = element_blank()) +   # Remove y-axis text (labels)
  labs(title = "DMBT1 Gene Expression Across Developmental Stages", 
       x = "Developmental Stage", 
       y = "DMBT1 Expression")

# Print the plot
print(DMBT1)

#G3bp2
normal<- read.csv("larvalnormalizedfilteredreads.csv", row.names=1)
normalfil<- subset(normal,rownames(normal)%in%c("TRINITY_DN38058_c0_g1_i1"))
normaltransform<- as.data.frame(t(normalfil))
normaltransform
rownames_to_column(normaltransform, var = "sequenceID")
normaltransform$sequenceID<- 38:84
mergedbby<- merge(normaltransform, meta, by="sequenceID")
mergedbby
# Filter out rows where developmental_stage is "planula_96h_F" or "planula_72h_F"
mergedbby <- mergedbby %>%
  filter(!developmental_stage %in% c("planula_96h_UF", "planula_72h_UF"))

mergedbby$developmental_stage <- factor(mergedbby$developmental_stage, 
                                        levels = c("zygote", "2_cell", "4_cell", "blastula", "gastrula", 
                                                   "planula_48h", "planula_72h_F","planula_96h_F"))
mergedbby

G3bp2 <- ggplot(data = mergedbby, aes(x = developmental_stage, y = TRINITY_DN38058_c0_g1_i1)) + 
  geom_smooth(aes(group = 1), method = "loess", se = FALSE) + 
  geom_point() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_blank(),    # Remove external ticks on the y-axis
        axis.text.y = element_blank()) +   # Remove y-axis text (labels)
  labs(title = "G3bp2 Gene Expression Across Developmental Stages", 
       x = "Developmental Stage", 
       y = "G3bp2 Expression")

# Print the plot
print(G3bp2)

ggsave("~/Desktop/LarvalDevRI22/Figures/G3bp2acrossstagesFED.pdf", G3bp2, height = 4, width = 8)

###GOMWU on differentially expressed immune genes 
##Parsing files first 
TopHits <- read.table("astrangia_annos_clean.txt", header = TRUE)
TopHits[1:10,1:10]
TopHits=TopHits[,c(1:2)]
TopHits=separate(data=TopHits, col=V2, into = c("sp","spID","gene_name", sep="\\|"))
TopHits=TopHits[,c(1,3)]


##Parsing DESeq Files to get input files for Matz Script
Zygote2cellCSV <- read.csv(file = "~/Desktop/LarvalDevRI22/Interceptdevelopmental_stage4_cell.csv")
Zygote2cellCSV[1:10,1:7]
ParsedZygote2cell=Zygote2cellCSV[,c(1,3)]
ParsedZygote2cell[1:10,1:2]
write.csv(ParsedZygote2cell, file = "Parsed_2cell4cell.csv", row.names = FALSE)


setwd("~/Desktop/Masters/Spring_2022/FuessLab/Transcriptomics")
if(!require(ape))
  install.packages("ape")
library(ape)

# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.


# Edit these to match your data file names: 
input="Parsed_2cell4cell.csv" # two columns of comma-separated values: Trinity id and log 2 fold change
goAnnotations="merged_GO_data.tab" # two-column, tab-delimited, one line per gene, multiple GO IDs separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="GO_hierarchy.obo" # make sure this is in the directory for it to work-- download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="/usr/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


##Output file will be dissim_CC_"name of input file".tab - I stopped here for my thesis!

# ----------- Plotting results

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  

# text representation of results, with actual adjusted p-values
results[[1]]


# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
# plot(results[[2]],cex=0.6)
# abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs

