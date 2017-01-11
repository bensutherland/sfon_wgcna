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

setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/sfon_wgcna")
setwd("/Users/wayne/Documents/bernatchez/Sfon_projects/SfeQ/eQTL_analysis/edgeR-norm")

# import count data from each sample file:
files <- list.files(path="/Users/wayne/Documents/bernatchez/Sfon_projects/SfeQ/eQTL_analysis/edgeR-norm/01_input_data/", 
                    pattern="*htseq_counts.txt.tab")

my.counts <- readDGE(files, path="/Users/wayne/Documents/bernatchez/Sfon_projects/SfeQ/eQTL_analysis/edgeR-norm/01_input_data")
  # note that there is an argument for DGEList 'genes' to provide annotation for tags/transcripts/genes
str(my.counts)
dim(my.counts) #total unique tags = 69440

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

