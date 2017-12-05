# Identify which transcripts belong to contiguous blocks and identify unique blocks

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Import data
data <- read.table(file = "sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bed")
data <- data[,c(1:4)]
colnames(data) <- c("target.contig","start","end","transcript")
str(data)


# add shared gene column to populate with counter
shared.gene <- rep(NA, times=nrow(data)) 
data <- cbind(data, shared.gene)
str(data)

# testing
#data <- data[1:200,]


### Set nulls
chr <- NULL; prev.chr <- NULL
start <- NULL ; prev.start <- NULL
end <- NULL ; prev.end <- NULL
output <- NULL

# Set up for false start
prev.chr <- 1
prev.end <- 0

# Set counter
counter <- 1

# Loop
for(i in 1:nrow(data)){
  
  # set the current variables for this round
  chr <- data$target.contig[i]
  start <- data$start[i]
  end <- data$end[i]
  
  # add redundant statement if same chr and start before prev end
  if(chr == prev.chr && start < prev.end ){
    data$shared.gene[i] <- counter
    print("redundant")
  } else {
    counter <- counter + 1
    data$shared.gene[i] <- counter
    print("new.gene")
    
    # if not redundant, set previous variables
    prev.chr <- chr
    prev.start <- start
    prev.end <- end
  }
}

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




# Import the wgcna object (female)
geneinfo <- read.delim2(file = "female_geneInfo_added_annot.txt", header = T)
colnames(geneinfo)
head(geneinfo)

#### Find best one to keep
head(data.lengths.sorted, n = 10)

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

write.table(x = result, file = "single_transcript_per_gene.txt", quote = F, col.names = T, row.names = F, sep = "\t")



### Get proportions of each





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
length(unique(x = data.lengths$shared.gene))

max(data.lengths$shared.gene)
# NULLS
