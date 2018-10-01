## Differential Expression Analysis (sex-biased expression)
# This is the step after normalizing, but before building networks.

# Inputs: The RData object from the script '01_edgeR_normalization'
# 1) save.image(file = "02_input_data/sfon_wgcna_01_output.RData") # save out existing data 
#rm(list=ls())

## Install Packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
require("edgeR")

#biocLite("locfit")
require("locfit")

# Set working directory
#setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/sfon_wgcna") #wayne
setwd("~/Documents/10_bernatchez/01_sfon_eqtl/sfon_wgcna") #logan

load("02_input_data/sfon_wgcna_01_output.RData")

##### To Do #####
# 1. Expressed / Not Expressed analysis
# 2. More appropriate low expression filter in terms of number of samples needed to retain transcript?


# The DGElist is here (note this is not cpm nor log2)
my.counts
# (Dispersions have already been estimated)


#### Remove unneeded samples for DE analysis ####
# Subset out the Arctic Charr
ac.pos <- grep(pattern = "SRR", x = rownames(my.counts$samples), perl = T)
my.counts <- my.counts[, -c(ac.pos), keep.lib.sizes=T]
my.counts$samples
dim(my.counts$samples)

# Subset out the parents
parent.pos <- grep(pattern = "101|102", x = rownames(my.counts$samples), perl = T)
parent.pos
my.counts <- my.counts[, -c(parent.pos), keep.lib.sizes=T]
dim(my.counts$samples)

#### Identify phenotypic sex for samples ####
## Get sex info per sample
# Get sample names of dgelist
dgelist.samples.df <- as.data.frame(rownames(my.counts$samples)) # get sample names
colnames(dgelist.samples.df) <- "file.name"
head(dgelist.samples.df)
dim(dgelist.samples.df) # 100 samples

# Combine w/ sex phenos
dgelist.ordered.interp <- merge(x = dgelist.samples.df, y = interp, by = "file.name", 
      sort = F) # this is essential, will keep in same order as the dgelist
# As check, these two should match
head(dgelist.ordered.interp$file.name)
head(rownames(my.counts$samples))

# Make binary vector with the sex of the samples in the order of the dgelist
levels(dgelist.ordered.interp$sex) <- c(0,1) # females 0, male 1
dgelist.ordered.interp$sex

# Add this binary sex vector as the 'group' for the DGElist
my.counts$samples$group <- dgelist.ordered.interp$sex # now the group is by sex

# Again, make sure these match (file.name and sex)
head(my.counts$samples, n = 10)
head(dgelist.ordered.interp[, c("file.name", "sex")], n = 10)


#### DE Analysis ####
design <- model.matrix(~my.counts$samples$group)
colnames(design)[2] <- "sex"
fit <- glmFit(y = my.counts, design = design)
# cont.matrix <- makeContrasts(sex, levels=design) # intercept is female
# ?contrasts.fit
lrt <- glmLRT(fit)
topTags(lrt, n=50)

result <- topTags(lrt, n = 1000000) # extract information from the glmLRT
sum(result$table$FDR < 0.05) # very similar to that approach below, slightly diff b/c glmLRT?
sum(result$table$FDR < 0.05 & abs(result$table$logFC) > log2(1.5))

# Export
result.output <- result$table # put into easily exported result
dim(result.output)
result.output$transcript <- rownames(result.output)
head(result.output)
result.output <- result.output[,c("transcript", "logFC", "logCPM", "LR", "PValue", "FDR")]
head(result.output)

# Limit by FC and FDR
result.output.filt <- result.output[abs(result.output$logFC) > log2(1.5) & result.output$FDR < 0.05 ,]
dim(result.output.filt)

# Export the results
write.csv(x = result.output.filt, file = "04_results/transcriptome_DE_results_0.05_fc1.5.csv", quote = F, row.names = F)


#### Annotation ####
annot <- read.table(file = "02_input_data/sfontinalis_contigs_annotation_report_v1.0_shortform.txt"
                    , sep = "\t", header = T)

head(annot)
head(result.output) # all transcript analyzed
background.annot <- merge(x = result.output, y = annot, by.x = "transcript", by.y = "transcript_id")
head(background.annot)

filt.annot <- merge(x = result.output.filt, y = annot, by.x = "transcript", by.y = "transcript_id")
head(filt.annot)


#### Result Exploration ####
high.fc.thresh.log2 <- log2(4)
head(filt.annot)
moderate.female.bias <- filt.annot[filt.annot$logFC < 0 & filt.annot$logFC > -c(high.fc.thresh.log2),]
high.female.bias <- filt.annot[filt.annot$logFC < -c(high.fc.thresh.log2),]
moderate.male.bias <- filt.annot[filt.annot$logFC > 0 & filt.annot$logFC < c(high.fc.thresh.log2),]
high.male.bias <- filt.annot[filt.annot$logFC > c(high.fc.thresh.log2),]

head(high.female.bias)

datasets <- c("moderate.male.bias", "high.male.bias", "moderate.female.bias", "high.female.bias")

# how many genes in each?
for(i in 1:length(datasets)){
  print(datasets[i])
  print(paste("number transcripts = ", nrow(get(datasets[i]))))
  
  # export result
  filename <- paste0("04_results/", datasets[i], "_DEGs.txt")
  write.csv(x = get(datasets[i]), file = filename, quote = F, row.names = F)
  
  # Give reports
  print(paste("median logCPM = ", median(get(datasets[i])$logCPM)))
  print(paste("mean logCPM = ", mean(get(datasets[i])$logCPM)))
}


# What about the general data
dim(result.output)
median(result.output$logCPM)
mean(result.output$logCPM)

# Need to now figure out if the low expression is leading to a problem in CPM

# It looks like the format is 1 / 0, where positive log2 FC is 1 biased (male biased) 

#### Plot single transcripts
# for example:
TOI <- "QSF_ABCC5.1.1"
boxplot(my.counts$counts[TOI,] ~ my.counts$samples$group)

# OTHER APPROACH
# from https://gist.github.com/jdblischak/11384914
et <- exactTest(my.counts)
results_edgeR <- topTags(et, n = 100000, sort.by = "none")
str(results_edgeR)
head(results_edgeR$table)
sum(results_edgeR$table$FDR < .05) # how many genes w/ FDR < cutoff
plotSmear(et, de.tags = rownames(my.counts$counts)[results_edgeR$table$FDR < 0.05])
abline(h = c(-2, 2), col = "blue")


# More data exploring and plotting
plotMDS(x = my.counts, cex= 0.8, labels = my.counts$samples$group)
