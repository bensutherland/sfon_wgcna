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




