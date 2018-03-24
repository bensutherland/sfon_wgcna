# Plots gene segments on chromosomes

# Inputs: 
## A. chr lengths: chr_of_interest_genome_lengths.txt
## B. txome to genome align file w modules: <sex>_single_transcript_per_gene.txt
## C. transcript/chr names of interest: <sex>_fisher_test_chr_enrich.txt

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Choose sex of coexpression modules
sex <- "female"
#sex <- "male"

##### 1. Input files ####
# A. chr length file
chr.length <- read.delim2(file = "chr_of_interest_genome_lengths.txt", header = F, sep = "\t")
head(chr.length)

# B. txome to genome align file w/ modules
in.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
data.file <- read.delim2(file = in.filename, header = T, sep = "\t")
head(data.file)
#issue w/ first line as NA

# C. transcripts/chr of interest file
in.filename <- paste(sex, "_fisher_test_chr_enrich.txt", sep = "")
chr.mods <- read.delim(file = in.filename, sep = "\t", header = T)
head(chr.mods)


#### 2. Subset to the chr/mod of interest #### 
pval.cutoff <- 0.01

print(c("You are analyzing dataset: ", sex, "with enrich pval cutoff of: ", pval.cutoff ), sep = "", quote = F)
sig.chr.mods <- chr.mods[chr.mods$pval < pval.cutoff, ] # only show those enriched with p < 0.01
sig.chr.mods

# Merge signif data w/ chr length
sig.chr.mods.and.length <- merge(x = sig.chr.mods, y = chr.length, by.x = "chr", by.y = "V1")
colnames(sig.chr.mods.and.length)[6] <- "length"  # issue make this more adaptive
head(sig.chr.mods.and.length)

# Remove low.corr as these won't have modules
sig.chr.mods.and.length <- sig.chr.mods.and.length[sig.chr.mods.and.length$mod != "low.corr", ]

# Remove un-enriched as these are not of interest (keep only those w/ at least 3 genes)
sig.chr.mods.and.length <- sig.chr.mods.and.length[sig.chr.mods.and.length$g.in.chr.in.mod > 2, ]

# Remove grey as this is not a true module
sig.chr.mods.and.length <- sig.chr.mods.and.length[sig.chr.mods.and.length$mod != "grey", ]

# sort it
# # Make the data ordered so that it matches the y-axis
sig.chr.mods.and.length <- sig.chr.mods.and.length[with(sig.chr.mods.and.length, order(sig.chr.mods.and.length$chr)), ]



#### 3. Plot ####
### Set up Plot ###
# Find the necessary ylim for the plot
max.ylim <- length(sig.chr.mods.and.length$chr) # based on the total unique chr involved

# Find the necessary xlim for the plot
max.xlim <- max(sig.chr.mods.and.length$length) # based on the longest chr

# Draw empty plot
par(mfrow=c(1,1), mar = c(4,7,3,6), cex = 0.6)

# Plot empty graph
plot(1, type="n", xlab="Position (bp)", ylab="", xlim=c(1, max.xlim), ylim=c(0, max.ylim), yaxt = "n")

# Use the chromosome names instead of the fasta accessions
lbls0 <- gsub(pattern = "NC_0273|\\.1", replacement = "", x =  sig.chr.mods.and.length$chr, perl = T)
lbls <- as.numeric(lbls0) + 1
lbls

axis(2, at=1:max.ylim
     , labels=paste("chr_", lbls, sep = "")
     #, labels=paste("chr_", seq(1:10))
     , las = 1)

axis(4, at=1:max.ylim
     , labels=sig.chr.mods.and.length$mod
     , las = 1)

### Loops to plot positions of genes in the chromosome
sig.chr.mods.and.length # this is our data to use

# Pairs we want to plot
plot.data <- unique(paste(sig.chr.mods.and.length$chr, sig.chr.mods.and.length$mod, sep = "-")) #join and uniq
plot.data

pairs.of.interest <- strsplit(plot.data, "-") # split it back up into constituents in a list
pairs.of.interest 

# Loop to plot positions of genes within the chromosome
chr.of.interest <- NULL; mod.of.interest <- NULL
sorted.sig <- NULL
chr.pos <- NULL; chr.max <- NULL

for(i in 1:length(pairs.of.interest)){
  # print(pairs.of.interest[[i]])
  chr.of.interest <- pairs.of.interest[[i]][1]
  mod.of.interest <- pairs.of.interest[[i]][2]
  print(chr.of.interest)
  print(mod.of.interest)
  
  # Subset data to find the genes and their positions for the chr and mod of interest
  data.of.interest <- subset(data.file, target.contig==chr.of.interest & moduleColor==mod.of.interest)
  print(data.of.interest)
  
  # # Make the data ordered so that it matches the y-axis
  # sorted.sig <- sig.chr.mods.and.length[with(sig.chr.mods.and.length, order(sig.chr.mods.and.length$chr)), ]
  
  # Find the max length of the chr of interest, and the y position for the y-axis
  #chr.pos <- which(sig.chr.mods.and.length$chr == chr.of.interest) ## ISSUE HERE
  chr.pos <- i
  
  chr.max <- sig.chr.mods.and.length[sig.chr.mods.and.length$chr == chr.of.interest, "length"][1]
  # note that the [1] is used as chr are repeated in the file that is being queried
  
  # inner Loop to add segments
  start <- NULL; stop <- NULL
  chr.max
  
  for(g in 1:length(data.of.interest$transcript)){
    start <- data.of.interest$start[g]
    stop <- data.of.interest$end[g]
    #print(c(start, stop))
    segments(x0 = start, y0 = chr.pos, x1 = stop, y1 = chr.pos, col = mod.of.interest, lwd = 3)
    points(x = start, y = chr.pos, pch="/")
    points(x = stop, y = chr.pos, pch="|")
    points(x = chr.max, y = chr.pos, pch=23)
  }
  
  
}


# save out as
# <sex>_chr_mod_enrichment.pdf in 5 x 3 inch


