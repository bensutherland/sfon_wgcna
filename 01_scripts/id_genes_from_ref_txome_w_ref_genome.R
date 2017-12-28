# Identify putative unique genes from a reference transcriptome
# by using a reference genome
# input: bed file; NOTE: MUST BE SORTED BY "target.contig" THEN "start"

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Set filenames
bed.filename <- "sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bed"
transcript.lengths.filename <- "sfontinalis_contigs_unwrap_seq_lengths.txt"
# Also need: male_geneInfo_added_annot.txt OR female_geneInfo_added_annot.txt 

#### 0. Import data ####
# Import bed file (ref txome against genome)
data <- read.table(file = bed.filename)
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

# Prevent error due to NULLs on initial run
prev.chr <- 1 ; prev.end <- 0

# Initialize counter
counter <- 1

# Identify which genes are in continuous blocks
for(i in 1:nrow(data)){
  
  ## Set the current variables for this round
  chr <- data$target.contig[i]
  start <- data$start[i]
  end <- data$end[i]
  
  # debugging
  #print(c(i, end))
  
  # add same ID if transcript aligns to same chr and starts before the prev. end
  if(chr == prev.chr && start < prev.end ){
    data$shared.gene[i] <- counter
    
    # Find maximum value in "end" column for this specific shared.gene ID
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

#### 2. Incorporate fasta record lengths ####

# Attach fasta record lengths to this file
# Bring in fasta record length file
transcript.lengths <- read.table(file = transcript.lengths.filename
                            , col.names = c("transcript", "length.bp"))
str(transcript.lengths)

# Merge lengths file with non-redundant transcript file
data.lengths <- merge(x = data, y = transcript.lengths, by = "transcript", sort = F)

str(data.lengths)
colnames(data.lengths)
head(data.lengths)

# Sort data.lengths by shared.gene ID then by length
data.lengths.sorted <- data.lengths[with(data.lengths, order(data.lengths$shared.gene, data.lengths$length.bp
                          , decreasing= T)), ]

# Add empty vector to be populated with whether the gene is expressed or not in the expression data (to be imported)
is.present <- rep("NA", times = length(data.lengths.sorted$transcript))
data.lengths.sorted.choose <- cbind(data.lengths.sorted, is.present )

#### 3. Import gene expression object ####
# Choose sex of coexpression modules
#sex <- "female"
sex <- "male"

# Create filename
geneInfo.filename <- paste(sex, "_geneInfo_added_annot.txt", sep = "")

# Import WGCNA object
geneinfo <- read.delim2(file = geneInfo.filename, header = T)
colnames(geneinfo)

# Determine which transcripts were expressed in the selected WGCNA object
# extract necessary data (transcript ID and module Color)
expressed.transcripts <- geneinfo[,c("transcript_id","moduleColor")]
head(expressed.transcripts)                                 
colnames(expressed.transcripts) <- c("transcript", "moduleColor")
                                     
# Merge with the rest of the data (keeping the data file for those that don't have an expr transcript)
all.data <- merge(x = data.lengths.sorted.choose, y = expressed.transcripts, by = "transcript", all.x = T)
dim(all.data)
head(all.data)

# Sort again, by shared.gene ID, then by length
all.data <- all.data[with(all.data, order(all.data$shared.gene, all.data$length.bp
                                                             , decreasing= T)), ]
head(all.data, n = 10)
str(all.data)

# Use the is.present column, first make it a character
all.data$is.present <- as.character(all.data$is.present)
# Then if there is a module color (not NA), give it 'yes' and if no module color (NA), give it 'z'
# The z is used for ease of alphabetical ordering
all.data[is.na(all.data$moduleColor)==F, "is.present"] <- "yes" 
all.data[is.na(all.data$moduleColor)==T, "is.present"] <- "no" 
head(all.data)

# Sort again, by gene, is.present, then length
all.data <- all.data[with(all.data, 
                          order(all.data$shared.gene, all.data$is.present, all.data$length.bp, decreasing = T)
                          ), ]
tail(all.data, n = 10)

#Testing
# all.data <- all.data[1:20,]
# head(all.data)
# dim(all.data)

#### 4. Obtain one record per shared.gene ID ####
# Set nulls
goi <- NULL; section <- NULL; result <- NULL

# Make a vector of unique gene names
unique.genes <- unique(all.data$shared.gene)

# 
for(i in unique.genes){
  # identify the name of this round's goi
  goi <- unique.genes[i]
  #print(goi)
  # select the rows that are within the goi
  section <- all.data[all.data$shared.gene==goi,]
  
  # take only the first one (because it's been sorted this works)
  result <- rbind(result, head(section, n = 1))
}

# This provides a dataset of all the top genes per goi.

head(result)
length(result$transcript)
length(unique(result$transcript)) # issue here

# Save out
out.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
write.table(x = result, file = out.filename, quote = F, col.names = T, row.names = F, sep = "\t")












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
