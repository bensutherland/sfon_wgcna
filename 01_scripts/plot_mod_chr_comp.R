# Produce graphs of the composition of each module in terms of chromosomes of containing genes

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Choose sex of coexpression modules
sex <- "female"
#sex <- "male"

# Import
in.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
results <- read.delim2(file = in.filename, header = T, sep = "\t")
head(results)

# remove row NA
results <- results[-c(which(is.na(results$transcript))),]

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

#### 01c. Baseline: genes per chromosome of interest ####
background.numbers.sorted[background.numbers.sorted$target %in% chromosomes.of.interest,]


#### 02. Count genes per chr for each module #### 

# Per module, count genes in each chromosome of interest
current.module.set <- NULL; result.list <- list() ; genes.per.chr.per.mod.list <- list()
info.set <- NULL ; info.set.all <- NULL

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
    #names(genes.per.chr.per.mod.list[[m]] <- chromosomes.of.interest)
  }
}

# summary data:
info.set.all # gives info on which modules contain how many transcripts
sum(as.numeric(info.set.all[c(-1,-3),2])) # to count up



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
par(mfrow=c(1,1), mar = c(0,4,2,0), cex = 0.6)

# Set up variables if composite plotting
plot.filename <- paste(sex, "_modules_by_chromosomes.pdf")
pdf(file = plot.filename, width= 12, height = 10)
# Importantly, set the parameters based on sex (number modules)
sex.par <- list() ; sex.par[["female"]] <- c(4,4); sex.par[["male"]] <- c(7,4)
par(mfrow=sex.par[[sex]], mar = c(0,4,2,0), cex = 0.6)


# Make a baseline pie chart (not sep by module)
slices <- c(background.numbers.sorted$num.genes[1:29])
lbls <- background.numbers.sorted$target[1:29]
lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = lbls, perl = T)
pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls,"-", pct, "%", sep="")

pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
    , main = paste("Baseline with", sum(slices), "genes"))

# Then plot separately per module
for(i in 2:length(genes.per.chr.per.mod.list)){
  slices <- genes.per.chr.per.mod.list[[i]]
  slices <- slices[slices!=0] # don't plot any slices that are empty (required with THIS)
  # lbls <- chromosomes.of.interest
  lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = chromosomes.of.interest, perl = T)
  pct <- round(slices/sum(slices)*100)
  lbls2 <- paste(lbls,"-", pct, "%", sep ="")
  lbls2 <- lbls2[grep(pattern = "-0%", x = lbls2, invert = T)] # don't plot any labels when empty (required with THIS)
  module.this.round <- names(genes.per.chr.per.mod.list)[i]
  
  pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
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
# Write out results
filename.output <- paste(sex, c("count_of_genes_per_chr_per_module.txt"), sep = "_")
write.table(x = all.output
            , file = filename.output
            , quote = F
            , sep = "\t"
            , row.names = F)





### Leftover Code
# How to subset based on module
# m <- "grey"
# head(results[results$moduleColor %in% m,])

# This counts up how many instances of each chromosome are within the current module
#aggregate(x = current.module.set$transcript, by = list(target = current.module.set$target.contig), FUN = length)
