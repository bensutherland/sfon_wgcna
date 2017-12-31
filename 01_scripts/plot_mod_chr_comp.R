# Produce graphs of the composition of each module in terms of chromosomes of containing genes

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Choose sex of coexpression modules
#sex <- "female"
sex <- "male"

number.modules.by.sex <- c(0,28)
names(number.modules.by.sex) <- c("female","male")

number.modules.by.sex[sex]

# Import
in.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
results <- read.delim2(file = in.filename, header = T, sep = "\t")
head(results)

# remove NA
results <- results[-c(which(is.na(results$transcript))),]

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

# Collect data
# background data
background.numbers.sorted[background.numbers.sorted$target %in% chromosomes.of.interest,]

# foreground data



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
sex.par <- list() ; sex.par[["female"]] <- NULL; sex.par[["male"]] <- c(7,4)
par(mfrow=sex.par[[sex]], mar = c(0,4,2,0), cex = 0.6)


# Make a baseline pie (not sep by module)
slices <- c(background.numbers.sorted$x[1:29])
lbls <- background.numbers.sorted$target[1:29]
lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = lbls, perl = T)
pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls,"-", pct, "%", sep="")

pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
    , main = paste("Baseline with", sum(slices), "genes"))

# Then plot separately per module
for(i in 2:length(test.list)){
  slices <- test.list[[i]]
  slices <- slices[slices!=0] # don't plot any slices that are empty (required with THIS)
  # lbls <- chromosomes.of.interest
  lbls <- gsub(pattern = "NC_0273|\\.1", replacement = "", x = chromosomes.of.interest, perl = T)
  pct <- round(slices/sum(slices)*100)
  lbls2 <- paste(lbls,"-", pct, "%", sep ="")
  lbls2 <- lbls2[grep(pattern = "-0%", x = lbls2, invert = T)] # don't plot any labels when empty (required with THIS)
  module.this.round <- names(test.list)[i]
  
  pie(x = slices, labels = lbls2, col = palette[1:length(lbls2)]
      , main = paste(module.this.round, "with", sum(slices), "genes"))
}

dev.off()


### Save out values




### Leftover Code
# How to subset based on module
# m <- "grey"
# head(results[results$moduleColor %in% m,])

# This counts up how many instances of each chromosome are within the current module
#aggregate(x = current.module.set$transcript, by = list(target = current.module.set$target.contig), FUN = length)
