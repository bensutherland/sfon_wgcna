# Plot genes on chromosomes
# input: chromosome lengths, bed file, transcript names of interest

# use: <sex>_single_transcript_per_gene.txt, has start, end and moduleColor, transcript and target contig

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")

# Choose sex of coexpression modules
#sex <- "female"
sex <- "male"

### Input files
# 
in.filename <- paste(sex, "_single_transcript_per_gene.txt", sep = "")
data.file <- read.delim2(file = in.filename, header = T, sep = "\t")

head(data.file)

# chr length file
chr.length <- read.delim2(file = "chr_of_interest_genome_lengths.txt", header = F, sep = "\t")
head(chr.length)


### HERE ###


par(mfrow=c(1,1), mar = c(4,4,3,3), cex = 0.6)

# use plot to draw empty graph
plot(1, type="n", xlab="", ylab="", xlim=c(1, 10), ylim=c(0, 10), yaxt = "n")
## note: replace xlim with the max length of chromosomes of interest
axis(2, at=1:10, labels=paste("chr_", seq(1:10)), las = 1)


# Plot multiple chromosomes this way:
segments(x0 = 1, y0 = 10, x1 = 10, y1 = 10, col = "black", lwd = 1)

# add genes this way
chr.pos <- 10 # this would be the level that the chr is at on the graph
start <- 2 # this would be the starting position of the gene
stop <- 3 # this would be the ending position of the gene

segments(x0 = start, y0 = chr.pos, x1 = stop, y1 = chr.pos, col = "blue", lwd = 5)


# Plot a couple examples..

# e.g. lightsteelblue1, contig NC_027317.1 
#or NC_027312.1	yellowgreen
#or NC_027302.1	grey60

#or NC_027302.1	yellowgreen
# NC_027312.1	skyblue
# NC_027324.1	turquoise
# NC_027325.1	yellow
# NC_027305.1	plum1
# NC_027304.1	steelblue

## these ones:
#chr	mod	pval	g.in.chr.in.mod	g.not.in.chr.in.mod
# NC_027317.1	lightsteelblue1	1.41506678012974E-05	4	2
# NC_027312.1	yellowgreen	0.0006309297	5	10
# NC_027316.1	turquoise	0.0018427814	1	348
# NC_027302.1	grey60	0.0021096624	4	9
# NC_027305.1	low.corr	0.003513191	372	9371
# NC_027302.1	yellowgreen	0.0037478207	4	11
# NC_027312.1	skyblue	0.0042669006	4	10
# NC_027324.1	turquoise	0.0050767741	16	333
# NC_027325.1	yellow	0.0052164004	5	55
# NC_027305.1	plum1	0.0061440002	4	14
# NC_027304.1	steelblue	0.0073398353	4	15
# NC_027305.1	grey	0.0087481258	348	6992
# NC_027320.1	sienna3	0.009282496	3	18

# choose those w/ p < 0.01, greater than 3 genes, and actual moduels (not low.corr or grey)
