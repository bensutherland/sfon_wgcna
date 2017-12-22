# This script identifies putative unique genes from a reference transcriptome
# that has been aligned against a reference genome to produce a bed file

#rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

#### 0. Import data ####
# Ref txome against genome bed file
data <- read.table(file = "sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bed")
colnames(data) <- c("target.contig", "start", "end", "transcript", "score", "sense")
head(data)
data <- data[,c(1:4)]

# Add empty vector to be populated with 'unique gene' counter
shared.gene <- rep(NA, times=nrow(data)) 
data <- cbind(data, shared.gene)
str(data)

#### 1. Link contiguous transcripts into genes ####

# testing
#data.bck <- data 
#data <- data.bck
#data <- data[1:200,]

### Set nulls
chr <- NULL; prev.chr <- NULL; start <- NULL ; prev.start <- NULL
end <- NULL ; prev.end <- NULL; output <- NULL

# Set up in case there is no initial values
prev.chr <- 1 ; prev.end <- 0

# Initialize counter
counter <- 1

# Loop
for(i in 1:nrow(data)){
  
  ## Set the current variables for this round
  chr <- data$target.contig[i]
  start <- data$start[i]
  end <- data$end[i]
  
  # debugging
  print(c(i, end))
  
  # add redundant indicator if transcript aligns to same chr, and starts before the previous end
  if(chr == prev.chr && start < prev.end ){
    data$shared.gene[i] <- counter
    
    # Find maximum end for this shared.gene
    prev.end <- max(data[which(data$shared.gene==counter), "end"])
    print(prev.end)
    
  } else {
    # increase the counter to signify a new 'shared.gene'
    counter <- counter + 1
    # and give this new 'shared.gene' value to this transcript
    data$shared.gene[i] <- counter
    # Then set new 'previous' variables
    prev.chr <- chr
    prev.start <- start
    prev.end <- end
  }
}


# Write out results
write.table(x = data
            , file = "sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted_w_redund_info.bed"
            , quote = F
            , sep = "\t"
            , row.names = F)

# View data
head(data)

# Attach fasta record lengths to this file
# Bring in fasta record length file
transcript.lengths <- read.table(file = "sfontinalis_contigs_unwrap_seq_lengths.txt"
                            , col.names = c("transcript", "length.bp"))
str(transcript.lengths)

# merge
data.lengths <- merge(x = data, y = fasta.lengths, by = "transcript", sort = F)

str(data.lengths)
colnames(data.lengths)
head(data.lengths)

# sort data.lengths
data.lengths.sorted <- data.lengths[with(data.lengths, order(data.lengths$shared.gene, data.lengths$length.bp
                          , decreasing= T)), ]

# Add column
is.present <- rep("NA", times = length(data.lengths.sorted$transcript))
data.lengths.sorted.choose <- cbind(data.lengths.sorted, is.present )



############# IMPORT GENE EXPRESSION OBJECT ###########
# Import the wgcna object (female)
geneinfo <- read.delim2(file = "female_geneInfo_added_annot.txt", header = T)
colnames(geneinfo)
head(geneinfo)

# Import the wgcna object (male)
geneinfo <- read.delim2(file = "male_geneInfo_added_annot.txt", header = T)
colnames(geneinfo)
head(geneinfo)

#### Find best one to keep
#head(data.lengths.sorted, n = 10)

# which transcripts were expressed?
colnames(geneinfo)
expressed.transcripts <- geneinfo[,c("transcript_id","moduleColor")]
head(expressed.transcripts)                                 
colnames(expressed.transcripts) <- c("transcript", "moduleColor")
                                     
all.data <- merge(x = data.lengths.sorted.choose, y = expressed.transcripts, by = "transcript", all.x = T)
dim(all.data)
head(all.data)

# sort again
all.data <- all.data[with(all.data, order(all.data$shared.gene, all.data$length.bp
                                                             , decreasing= T)), ]
head(all.data, n = 10)
str(all.data)


# Add a generic is.present column to the all.data file to sort upon without sorting on alphabetized colors in moduleColors
all.data$is.present <- as.character(all.data$is.present) # make the 'is.present' variable a character
str(all.data)

all.data[is.na(all.data$moduleColor)==F, "is.present"] <- "yes" # if present, say yes
all.data[is.na(all.data$moduleColor)==T, "is.present"] <- "z" # if present, say z
head(all.data)

# Sort again, this time by gene and is.present
all.data <- all.data[with(all.data, 
                          order(all.data$shared.gene, all.data$is.present, all.data$length.bp, all.data$length.bp)
                          ), ]
# this should make it so that the genes with moduleColor are selected first 
head(all.data, n = 10)
dim(all.data)


# Make a smaller version for testing
# all.data <- all.data[1:20,]
# head(all.data)
# dim(all.data)


#### Obtain a single record per 'gene' ####
# Set nulls
goi <- NULL; section <- NULL; result <- NULL

# make a vector of unique gene names
unique.genes <- unique(all.data$shared.gene)

for(i in unique.genes){
  # identify the name of this round's goi
  goi <- unique.genes[i]
  #print(goi)
  # select the rows that are within the goi
  section <- all.data[all.data$shared.gene==goi,]
  
  # take only the first one (because it's been sorted this works for the best)
  result <- rbind(result, head(section, n = 1))
}

# This provides a dataset of all the top genes per goi.

head(result)
dim(result)
length(unique.genes)

# write.table(x = result, file = "single_transcript_per_gene.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(x = result, file = "single_transcript_per_gene_male.txt", quote = F, col.names = T, row.names = F, sep = "\t")


##### Pie Charts #####
# Obtain number of transcripts per chromosome per module
#results <- read.delim2(file = "single_transcript_per_gene.txt", header = T, sep = "\t")
results <- read.delim2(file = "single_transcript_per_gene_male.txt", header = T, sep = "\t")
head(results)

modules <- unique(result$moduleColor)


# how to subset
m <- "grey"
head(results[results$moduleColor %in% m,])

# what chromosomes are the main ones?
#results$

# This counts up how many instances of each chromosome within the current module
#aggregate(x = current.module.set$transcript, by = list(target = current.module.set$target.contig), FUN = length)

# First characterize the number across all of the genes
background.numbers<- aggregate(x = results$transcript, by = list(target = results$target.contig), FUN = length)

# sort
background.numbers.sorted <- background.numbers[with(background.numbers, order(background.numbers$x, decreasing=T)), ]
head(background.numbers.sorted, n = 30)
# so the first 29 of those are the chromosomes...

chromosomes.of.interest <- background.numbers.sorted$target[1:29]

# make a pie ! (baseline)
slices <- c(background.numbers.sorted$x[1:29])
lbls <- background.numbers.sorted$target[1:29]
pct <- round(slices/sum(slices)*100)
lbls2 <- paste(lbls,"%", pct)

pie(x = slices, labels = lbls2, col = rainbow(length(lbls2)), main = "Baseline")


# What is the number of transcripts within module that are the same chromosome
current.module.set <- NULL; result.list <- list() ; test.list <- list()
#
for(m in modules){
  # subset dataset for each module color
  current.module.set  <- results[results$moduleColor %in% m,]
  print(dim(current.module.set))
  
  # now check how many copies of each 'chromosome.of.interest'
  for(c in chromosomes.of.interest){
    this.chr.count.this.mod <- length(grep(pattern = c, x = current.module.set$target.contig))
    test.list[[m]] <- c(test.list[[m]], this.chr.count.this.mod)
    #names(test.list[[m]] <- chromosomes.of.interest)
    }
}



### Next will be to produce charts of all of those..
for(i in 2:length(test.list)){
  slices <- test.list[[i]]
  lbls <- chromosomes.of.interest
  pct <- round(slices/sum(slices)*100)
  lbls2 <- paste(lbls,"%", pct, sep = "")
  module.this.round <- names(test.list)[i]
  
  pie(x = slices, labels = lbls2, col = rainbow(length(lbls2))
      , main = paste(module.this.round, "with", sum(slices), "genes"))
  # also add to the title the number of genes in the module
}


str(result.list)
result.list[["darkred"]]






##### HERE BE DRAGONS ###


# Get rid of wierd NA
# all.data$moduleColor <- gsub(all.data$moduleColor, pattern = "<NA>", replacement = "NA")
# Loop
# for(i in 1:nrow(data)){
#   
#   # set the current variables for this round
#   chr <- data$target.contig[i]
#   start <- data$start[i]
#   end <- data$end[i]
#   
#   # add redundant statement if same chr and start before prev end
#   if(chr == prev.chr && start < prev.end ){
#     data$shared.gene[i] <- counter
#     print("redundant")
#   } else {
#     counter <- counter + 1
#     data$shared.gene[i] <- counter
#     print("new.gene")
    
# First subset a bit...
#all.data.bck <- all.data
#all.data <- all.data[1:50,]





# 
# 
# 
# # Loop to 
# for(i in 1:nrow(all.data)){
#   #print(all.data$transcript[i])
#   
#   # if its an NA, don't do anything
#   if(is.na(all.data$moduleColor[i])==FALSE){
#     all.data$is.present[i] <- "yes"
#     #print("do nothing")
#   } 
#   # else { 
#   #   # if its not an NA, add a present flag
#   #   all.data$is.present[i] <- "yes"
#   # }
# }
# 
# # then add something that checks if the previous one of the same unit is taken..
#     
    

# if(chr == prev.chr && start < prev.end ){
#   data$shared.gene[i] <- counter
#   print("redundant")
# } else {
#   counter <- counter + 1
#   data$shared.gene[i] <- counter
#   print("new.gene")
# 
# 
# 


### FRAGMENTS
# how many unique gene designations?
# length(unique(x = data.lengths$shared.gene))
# 
# max(data.lengths$shared.gene)
# # NULLS
