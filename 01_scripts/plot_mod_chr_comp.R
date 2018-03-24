# Produce graphs of the composition of each module in terms of chromosomes of containing genes

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Choose sex of coexpression modules
#sex <- "female"
sex <- "male"

# Choose whether all modules or significant (sig) are to be plotted
select.type <- "all"
#select.type <- "sig"


# Import datafile of redundancy removed transcripts (sex-specific)
in.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
results <- read.delim2(file = in.filename, header = T, sep = "\t")
head(results)

# Remove the NA row
results <- results[-c(which(is.na(results$transcript))),]

# Convert moduleColor to character
results$moduleColor <- as.character(results$moduleColor)

#### 01a. Set up: module names ####

# what are the module names?
unique(results$moduleColor)

# rename moduleColor <NA> to low.corr, as these were not included in the WGCNA approach
results[which(is.na(results$moduleColor)), "moduleColor"] <- "low.corr"

# Define module names
modules <- unique(results$moduleColor)
modules

#### 01b. Set up: determine main chromosomes #####
# Count genes per genomic scaffold
background.numbers <- aggregate(x = results$transcript, by = list(target = results$target.contig), FUN = length)
names(background.numbers) <- c("target", "num.genes")

# Sort by largest containing chromosomes
background.numbers.sorted <- background.numbers[with(background.numbers, order(background.numbers$num.genes, decreasing=T)), ]
head(background.numbers.sorted, n = 30)
# Clearly, the first 29 of those are the chromosomes...

# Determine the main chromosomes of the genome
chromosomes.of.interest <- background.numbers.sorted$target[1:29]

# Write out chromosomes.of.interest file to use later
write.table(x = chromosomes.of.interest, file = "chr_of_interest_list.txt"
            , quote = F, sep = "\t"
            , row.names = F
            , col.names = F
            )
write.csv(chromosomes.of.interest, file = "chr_of_interest_list.csv", col.names = F)

#### 01c. Baseline: genes per chromosome of interest ####
# show the genes from all modules, low correlated and grey included, that were included in this sex-specific set of transcripts
background.numbers.sorted[background.numbers.sorted$target %in% chromosomes.of.interest,]

baseline.info <- background.numbers.sorted[background.numbers.sorted$target %in% chromosomes.of.interest,]
sum(baseline.info[,2]) # provide the total number of transcript locations in the chromosomes of interest

#### 02. Count genes per chr for each module #### 

# Per module, count genes in each chromosome of interest
current.module.set <- NULL; result.list <- list() ; genes.per.chr.per.mod.list <- list()
info.set <- NULL ; info.set.all <- NULL ; chr.order  <- NULL

for(m in modules){
  # Subset results per module, write to object 'current.module.set'
  current.module.set  <- results[results$moduleColor %in% m, ]
  #print(c("number transcript in this module is" , nrow(current.module.set)))
  
  # Collect information on the number of genes in this module and module name to use later
  info.set <- cbind(m, nrow(current.module.set)) # combine the name of the module with the number of transcripts
  info.set.all <- rbind(info.set.all, info.set)
  
  # Within module, how many genes per 'chromosome.of.interest'
  for(c in chromosomes.of.interest){
    this.chr.count.this.mod <- length(grep(pattern = c, x = current.module.set$target.contig))
    genes.per.chr.per.mod.list[[m]] <- c(genes.per.chr.per.mod.list[[m]], this.chr.count.this.mod)
  }
}

chromosomes.of.interest


# summary data:
info.set.all # gives info on which modules contain how many transcripts
sum(as.numeric(info.set.all[c(-1,-3),2])) # to count up


# Set output variable
expt <- paste(sex, select.type, sep = "_")


#### 03. Plot proportions of chr per mod ####
# Set up plotting
#install.packages("RColorBrewer")
library("RColorBrewer")

cols1 <- brewer.pal(n = 9, name = "Set1")
cols2 <- brewer.pal(n = 8, name = "Dark2")
cols3 <- brewer.pal(n = 10, name = "Set3")
cols4 <- brewer.pal(n = 11, name = "Spectral")
palette <- c(cols1,cols2,cols3,cols4)

# Plot
# If single plotting
# par(mfrow=c(1,1), mar = c(0,4,2,0), cex = 0.6)

# Set up variables if composite plotting
plot.filename <- paste(expt, "_modules_by_chromosomes.pdf")

# Set up size of output based on whether all or sig is to be shown
width.pdf <- list()
width.pdf[["female_all"]] <- 12 ; width.pdf[["male_all"]] <- 12
width.pdf[["female_sig"]] <- 10 ; width.pdf[["male_sig"]] <- 8

height.pdf <- list()
height.pdf[["female_all"]] <- 10 ; height.pdf[["male_all"]] <- 10
height.pdf[["female_sig"]] <- 3 ; height.pdf[["male_sig"]] <- 5


pdf(file = plot.filename, width= width.pdf[[expt]], height = height.pdf[[expt]])

# Importantly, set the parameters based on sex (number modules)

sex.par <- list()
sex.par[["female_all"]] <- c(4,5); sex.par[["male_all"]] <- c(7,4)
sex.par[["female_sig"]] <- c(1,3); sex.par[["male_sig"]] <- c(3,4)
par(mfrow=sex.par[[expt]], mar = c(0,4,2,2), cex = 0.6)


# Make a baseline pie chart (not sep by module)
slices <- c(background.numbers.sorted$num.genes[1:29])
lbls <- background.numbers.sorted$target[1:29]
lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = lbls, perl = T) # shorten name
# fix numbering issue (add 1 to each chromosome due to the difference in chr ID and accession ID)
lbls <- as.numeric(lbls) + 1

pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls,"-", pct, "%", sep="")

pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
    , main = paste("Baseline with", sum(slices), "genes"))

# Then plot separately per module
# set nulls
colors.subset <- NULL ; slices.plot <- NULL; lbls2.plot <- NULL; colors.subset.plot <- NULL


#### Get info on significant enrichment (new) ####
# Input transcripts/chr of interest file
in.filename <- paste(sex, "_fisher_test_chr_enrich.txt", sep = "")
chr.mods <- read.delim(file = in.filename, sep = "\t", header = T)
head(chr.mods)

# Only keep those that are significant
sig.chr.mods <- subset(chr.mods, pval < 0.01)
mods.to.plot <- as.character(unique(sig.chr.mods$mod))
# don't include low.corr or grey as these are not true modules
mods.to.plot <- mods.to.plot[mods.to.plot != "low.corr"]
mods.to.plot <- mods.to.plot[mods.to.plot != "grey"]

str(mods.to.plot)

# Build a loop that extracts the necessary elements from the larger list 'genes.per.chr.per.mod.list'
keep.info.list <- list()
mod <- c()
for(m in 1:length(mods.to.plot)){
  mod <- mods.to.plot[m]
  print(mod)
  keep.info.list[[mod]] <- genes.per.chr.per.mod.list[[mod]]
}

keep.info.list

### use an if statement to determine whether ALL or SELECT modules are to be plotted
if(select.type == "sig"){
  genes.per.chr.per.mod.list <- keep.info.list

} else{
  # don't change the main list
  print("Not changing the main list")
}



##### end new section ####

##### Plot #####

for(i in 1:length(genes.per.chr.per.mod.list)){
  slices <- genes.per.chr.per.mod.list[[i]]

  # lbls <- chromosomes.of.interest
  # Work with full data
  lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = chromosomes.of.interest, perl = T)
  
  # fix numbering issue (add 1 to each chromosome due to the difference in chr ID and accession ID)
  lbls <- as.numeric(lbls) + 1
  
  
  pct <- round(slices/sum(slices)*100)
  
  lbls2 <- paste(lbls,"-", pct, "%", sep ="")
  colors.subset <- palette[1:length(lbls2)]
  # These are all matched, and include zero value slices. need: slices, lbls2, colors.subset
  
  # Now remove slices and assoc. params (slice value, labels and colors)
  slices.plot <- slices[slices!=0] # don't plot any slices that are empty
  lbls2.plot <- lbls2[slices!=0]
  colors.subset.plot <- colors.subset[slices!=0]
  #lbls2 <- lbls2[grep(pattern = "-0%", x = lbls2, invert = T)] # don't plot any labels when empty
  
  
  # Set up plot
  module.this.round <- names(genes.per.chr.per.mod.list)[i]
  
  pie(x = slices.plot, labels = lbls2.plot, col = colors.subset.plot
      , main = paste(module.this.round, "with", sum(slices), "genes"))
}

dev.off()


### Save out values
# extract out gene counts from the list
genes.per.chr.per.mod.list # this contains all info
chromosomes.of.interest # this is the order of chromosomes

output <- NULL; output.temp <- NULL
for(i in 1:length(genes.per.chr.per.mod.list)){
  output.temp <- genes.per.chr.per.mod.list[[i]]
  output <- cbind(output, output.temp)
}
colnames(output) <- names(genes.per.chr.per.mod.list)
rownames(output) <- chromosomes.of.interest

head(output)

# add a row sum column for baseline
all.output <- cbind(output, rowSums(output))
colnames(all.output)[ncol(all.output)] <- "baseline"

all.output


#### 04. Save values data
#Write out results
filename.output <- paste(expt, c("count_of_genes_per_chr_per_module.txt"), sep = "_")
write.table(x = all.output
            , file = filename.output
            , quote = F
            , sep = "\t"
            , row.names = T
            , col.names = NA)





### Leftover Code
# How to subset based on module
# m <- "grey"
# head(results[results$moduleColor %in% m,])

# This counts up how many instances of each chromosome are within the current module
#aggregate(x = current.module.set$transcript, by = list(target = current.module.set$target.contig), FUN = length)
