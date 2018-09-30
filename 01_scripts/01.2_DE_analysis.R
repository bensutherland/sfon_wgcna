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
?contrasts.fit
lrt <- glmLRT(fit)
topTags(lrt, n=50)

result <- topTags(lrt, n = 1000000)
str(result)
sum(result$table$FDR < 0.05) # very similar to that approach below, slightly diff b/c glmLRT?
# in any case, this one can then be exported and used as DE analysis
write.table(result, file="04_results/transcriptome_DE_results.txt", sep = "\t")

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