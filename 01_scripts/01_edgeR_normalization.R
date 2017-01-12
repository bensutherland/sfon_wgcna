## Normalization of count data from HT-seq output
# This is the first step of the WGCNA repo, but the input data comes from:
# https://github.com/bensutherland/SE-reads_assemble-to-counts.git

# Install Packages
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
require("edgeR")

# User guides are here:
# browseVignettes("edgeR")
# edgeRUsersGuide()

# Set working directory
setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/sfon_wgcna")
# setwd("/Users/wayne/Documents/bernatchez/Sfon_projects/SfeQ/eQTL_analysis/edgeR-norm")

#### 1. Import Data ####
# Collect file names
files <- list.files(path="02_input_data/", pattern="*htseq_counts.txt")

# Import data
my.counts <- readDGE(files, path="02_input_data/", header = F) # Note: is an arg. for DGEList genes to provide annot.
str(my.counts)
dim(my.counts) #total unique tags = 69440

# Clean up file names
to.trim <- "_R1_trimmed.fastq.gz.bam_htseq_counts"
rownames(my.counts[[1]]) <- gsub(x = rownames(my.counts[[1]]), pattern = to.trim, replacement = "")
my.counts$samples$files <- gsub(x = my.counts$samples$files, pattern = to.trim, replacement = "")
dimnames(my.counts$counts)[[2]] <- gsub(x = dimnames(my.counts$counts)[[2]], pattern = to.trim, replacement = "")
str(my.counts) # much cleaner now

#### 2. Filter Data ####
#keep <- rowSums(my.counts[1,1] > 10) >= 2
# filter low expr tags
keep <- rowSums(cpm(my.counts)>0.5) >= 2
my.counts <- my.counts[keep,, keep.lib.sizes=FALSE] #keep.lib.sizes option is used to recom- pute the library sizes from the remaining tags
dim(my.counts) #now only have 1256 genes! #note: use 50cpm then get 2632 genes

# Normalization
my.counts <- calcNormFactors(my.counts, method = c("TMM"))
my.counts$samples #wow, note the norm factor on lib05!

# estimate dispersions (measure of inter-library variation for that tag)
my.counts <- estimateDisp(my.counts) # opts that did not work: , trend="none", robust=TRUE
summary(my.counts$prior.df) # this estimates the overall variability across the genome for this dataset
sqrt(my.counts$common.disp) #this gives the coeff of var of biological variation
plotBCV(my.counts)

# generate CPM matrix? #but does this use the normalized data i.e. calcNormFactors
test <- cpm(my.counts, normalized.lib.sizes=T, log=F)
head(test)
str(test)
write.csv(test, file = "test.csv")


#pca ?
plotMDS(my.counts)

