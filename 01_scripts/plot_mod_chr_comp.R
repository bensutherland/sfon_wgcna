# Produce graphs of the composition of each module in terms of chromosomes of containing genes

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Choose sex of coexpression modules
#sex <- "female"
sex <- "male"

# Import
in.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
results <- read.delim2(file = in.filename, header = T, sep = "\t")
head(results)

#### 01. Input data and determine genes per chromosome ####

# Define module names
modules <- unique(results$moduleColor)

# What chromosomes are the main ones?
# Find how many genes are found in each scaffold/chromosome
background.numbers <- aggregate(x = results$transcript, by = list(target = results$target.contig), FUN = length)

# Sort by largest containing chromosomes
background.numbers.sorted <- background.numbers[with(background.numbers, order(background.numbers$x, decreasing=T)), ]
head(background.numbers.sorted, n = 30)
# Clearly, the first 29 of those are the chromosomes...

chromosomes.of.interest <- background.numbers.sorted$target[1:29]


#### 02. Determine proportions of chromosomes per module #### 

# Determine the number of genes in each chromosome in each module
current.module.set <- NULL; result.list <- list() ; test.list <- list()
info.set <- NULL ; info.set.all <- NULL

for(m in modules){
  # subset dataset for each module color
  current.module.set  <- results[results$moduleColor %in% m,]
  #print(c("number transcript in this module is" , nrow(current.module.set)))
  
  info.set <- cbind(m, nrow(current.module.set)) # combine the name of the module with the number of transcripts
  info.set.all <- rbind(info.set.all, info.set) # prepare this for output later
  
  # Within this module, check how many copies of each 'chromosome.of.interest'
  for(c in chromosomes.of.interest){
    this.chr.count.this.mod <- length(grep(pattern = c, x = current.module.set$target.contig)) # for each c (chromosome)
    test.list[[m]] <- c(test.list[[m]], this.chr.count.this.mod)
    #names(test.list[[m]] <- chromosomes.of.interest)
  }
}

info.set.all # gives info on which modules contain how many transcripts
sum(as.numeric(info.set.all[c(-1,-3),2])) # to count up


#### 03. Plot proportions of chr per mod ####
# Set up plotting
#install.packages("RColorBrewer")
library("RColorBrewer")

cols1 <- brewer.pal(n = 9, name = "Set1")
cols2 <- brewer.pal(n = 8, name = "Set2")
cols3 <- brewer.pal(n = 10, name = "Set3")
cols4 <- brewer.pal(n = 11, name = "Spectral")
palette <- c(cols1,cols2,cols3,cols4)

# Plot
pdf(file = "modules_by_chromosomes.pdf", width= 12, height = 10)
par(mfrow=c(5,6), mar = c(0,4,2,0), cex = 0.6)

# Make a baseline pie (not sep by module)
slices <- c(background.numbers.sorted$x[1:29])
lbls <- background.numbers.sorted$target[1:29]
lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = lbls, perl = T)
pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls,"-", pct, "%", sep="")

pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)], main = "Baseline")

# Then plot separately per module
for(i in 2:length(test.list)){
  slices <- test.list[[i]]
  # lbls <- chromosomes.of.interest
  lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = chromosomes.of.interest, perl = T)
  pct <- round(slices/sum(slices)*100)
  lbls2 <- paste(lbls,"-", pct, "%", sep ="")
  module.this.round <- names(test.list)[i]
  
  pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
      , main = paste(module.this.round, "with", sum(slices), "genes"))
}

dev.off()


str(result.list)
result.list[["darkred"]]


pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
    , main = paste(module.this.round, "with", sum(slices), "genes"))

pie(x = slices[slices!=0]
    #, labels = lbls2
    , labels = lbls2[grep(pattern = "-0%", x = lbls2, invert = T)]
    , col = palette[1:length(lbls2)]
    , main = paste(module.this.round, "with", sum(slices), "genes"))


lbls2[grep(pattern = "-0%", x = lbls2, invert = T)]
?grep

#c hromosomes.of.interest <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = chromosomes.of.interest, perl = T)



### SCRAPS
# How to subset based on module
# m <- "grey"
# head(results[results$moduleColor %in% m,])

# This counts up how many instances of each chromosome are within the current module
#aggregate(x = current.module.set$transcript, by = list(target = current.module.set$target.contig), FUN = length)

# # Make a baseline pie (not sep by module)
# slices <- c(background.numbers.sorted$x[1:29])
# lbls <- background.numbers.sorted$target[1:29]
# pct <- round(slices/sum(slices)*100)
# lbls2 <- paste(lbls,"%", pct)
# 
# pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)], main = "Baseline")
# 
