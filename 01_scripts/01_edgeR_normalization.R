## Normalization of count data from eXpress output
# This is the first step of the WGCNA repo, but the input data comes from:
# https://github.com/bensutherland/Simple_reads_to_counts.git

# Inputs:
# 1) interpretation file
# 2) counts file (linear)

#rm(list=ls())

## Install Packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
require("edgeR")

biocLite("locfit")
require("locfit")

# User guides are here:
# browseVignettes("edgeR")
# edgeRUsersGuide()

# Set working directory
#setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/sfon_wgcna") #wayne
setwd("~/Documents/10_bernatchez/01_sfon_eqtl/sfon_wgcna") #logan

#### 1. Import Data, Set up Data ####
# Import interpretation file
interp <- as.data.frame(read.csv("00_archive/sfeq_interpretation_v1.csv"))
head(interp)
names(interp)


# Import counts file (i.e. the final output from Simple_reads_to_counts repo)
my.counts <- read.csv(file = "02_input_data/out.matrix.csv")


# Match the counts and interpretation
# fix my.counts colnames
new.names <- colnames(my.counts)
new.names <- gsub(new.names,pattern = ".eff.counts", replacement = "")
colnames(my.counts) <- new.names
names(my.counts)
dim(my.counts) #total unique tags = 69441
# my.counts.bk <- my.counts # just in case

# Fix filename in interp file
interp$file.name <- gsub(interp$file.name, pattern = "_trimmed.fastq.gz.bam_htseq_counts.trim.txt", replacement = "")
# note that we still don't have Arctic Charr in interpretation

# Set up DGEList
rownames(my.counts) <- my.counts[,1]
head(my.counts)

# Round all expression levels to no decimal place, remove transcript.id column (now rowname)
my.counts.round <- round(my.counts[,-1])
str(my.counts.round)


# create DGElist
my.counts <- DGEList(counts = my.counts.round)


#### 2. Filter Data ####
# Find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.reads.mapping.per.transcript <- 10
cpm.filt <- min.reads.mapping.per.transcript / min(my.counts$samples$lib.size) * 1000000
cpm.filt # min cpm filt

min.ind <- 5 # choose the minimum number of individuals that need to pass the threshold

# identify tags passing filter
keep <- rowSums(cpm(my.counts)>cpm.filt) >= min.ind # Find which transcripts pass the filter
table(keep) # gives number passing, number failing

# subset DGEList
my.counts <- my.counts[keep, , keep.lib.sizes=FALSE] #keep.lib.sizes = T retains original lib sizes, otherwise recomputes w remaining tags
dim(my.counts)


#### 3. Normalization ####
# Use TMM normalization, as it takes into account highly expressed genes that may take up sequencing rxn and make other genes look down-reg.
my.counts <- calcNormFactors(my.counts, method = c("TMM"))
my.counts$samples
plot(my.counts$samples$norm.factors ~ my.counts$samples$lib.size)

# #NOT WORKING AFTER CHANGED INPUT 2017-04-05
# my.counts$samples$files[my.counts$samples$norm.factors < 0.8] # find which files are outliers in terms of norm.factors
# # is there anything specific about these files?

# Estimate dispersions (measure inter-library variation per tag)
my.counts <- estimateDisp(my.counts) # note that this can use a design matrix when provided 
summary(my.counts$prior.df) # est. overall var. across genome for dataset
sqrt(my.counts$common.disp) #coeff of var, for biol. var
plotBCV(my.counts)


#### 4. Prepare Output ####
# generate CPM matrix
normalized.output <- cpm(my.counts, normalized.lib.sizes = TRUE, log= F)

# Compare the raw counts to the normalized cpm values (not log)
my.counts$counts[1:5, 1:5] # not normalized, raw counts
normalized.output[1:5, 1:5] # normalized lib size calculated cpm values

# output as normalized linear
write.csv(normalized.output, file = "03_normalized_data/normalized_output_matrix.csv")

# # output as normalized log2 (in progress)
normalized.output.log2 <- cpm(my.counts, normalized.lib.sizes = TRUE, log= T, prior.count = 1)
write.csv(normalized.output, file = "03_normalized_data/normalized_output_matrix_log2.csv")

# output object
save.image(file = "02_input_data/sfon_wgcna_01_output.RData") # save out existing data 


#### 5. Visualize data ####
#plot using sample IDs
plotMDS(x = my.counts, cex= 0.8) # note that this is supposed to be run on whatever you wrote calcNormFact() to

# #plot using sex
# plotMDS(x = my.counts, cex= 0.8
#         , labels = interp$sex[match(my.counts$samples$files, interp$file.name)])
# #plot using maturity and sex
# plotMDS(x = my.counts, cex= 0.8
#         , labels = paste(
#           interp$sex[match(my.counts$samples$files, interp$file.name)]
#             , interp$poids.sachet.foie[match(my.counts$samples$files, interp$file.name)]
#           , sep = ""))
# 
# # note, this is how matching works:
# interp$sex[match(my.counts$samples$files, interp$file.name)] # matches order 
# interp$sex #see not the same