# Plot sample dendrogram w/ traits, all polished

## Clean space
#rm(list=ls())

# Load setup wgcna from previous
load(file = "02_input_data/sfon_wgcna_setup.RData")

# Set up
datExpr0 <- datExpr0.bck

##### Remove AC samples #####
datExpr0 <- datExpr0[1:104, ] # remove the AC samples

#### Filter subset again specific to subset (from 03_WGCNA_BC_fem_mods_comp.R) ####
## Define which genes have > req num expressing samples
num.indiv <- 5

# Observe structure of data
dim(datExpr0)
datExpr0[1:5,1:5] 

# Create nulls
expressed <- NULL
keep.genes <- NULL
expressed.tally <- NULL

# loop over each gene (colnames)
for(i in 1:length(colnames(datExpr0))) { 
  
  # Note, here log2(1.19) is used as a proxy for cpm > 1.19 (from 01_edgeR_normalization.R)
  threshold <- log2(1.19)
  
  # How many samples have greater than threshold, per gene:
  expressed <- length(which(datExpr0[,i] > threshold))
  print(expressed)
  expressed.tally <- c(expressed.tally, expressed)
  
  # If the number of individuals for this gene is more than num.indiv above, put transcript in 'keep genes'
  if (expressed > num.indiv) {
    keep.genes <- c(keep.genes, colnames(datExpr0)[i])
  }
}

head(keep.genes)
length(keep.genes)


# # Give quantitative summaries of how many genes were expressed in how many samples
# length(rownames(datExpr0)) # number of samples
# 0.90 * length(rownames(datExpr0)) # 90% of samples
# table(expressed.tally  >  (0.90 * length(rownames(datExpr0)))) # provides answer
# 
# # Plot the freq of genes for each bin of number samples expressed?
# # filename <- paste("04_results/", REF, "_transcript_expr_in_num_samp.pdf", sep = "")
# # pdf(file = filename, width = 5, height = 4)
# 
# par(mfrow=c(1,1))
# hist(expressed.tally
#      , main = "Transcripts expressed in number samples"
#      , xlab = "Num. indiv. expr gene"
#      , ylab = "Num. transcripts in bin"
#      , las = 1) 
# 
# #dev.off()


# Filter by low expression
datExpr0.filt <- datExpr0[ ,keep.genes]
dim(datExpr0.filt)

# Replace orignal object with filtered
datExpr0 <- datExpr0.filt


#### 3. Data Quality Control ####
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #see tutorial if not true

# Plot samples to ID outliers (or remove Arctic Charr)
par(mfrow=c(1,1), mar = c(0,4,2,0), cex = 0.6)
sampleTree <- hclust(dist(datExpr0), method = "average") # default dist metric = euclidean; hclust agglomeration is "average"
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Remove AC samples
cutline <- 500

# To improve resolution on trends unrelated to the outlier cluster (i.e. liver weight), cut the tree.
abline(h = cutline, col = "red") # add to plot


# Identify which samples are in the largest portion of the cut tree
clust = cutreeStatic(sampleTree, cutHeight = cutline, minSize = 10)
table(clust) # clust 1 contains the samples we want to keep
keepSamples = (clust==1)

datExpr = datExpr0[keepSamples, ] # keep only main group
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #should be 104, note this includes parents 2x each


#### 4.a. Incorporate trait data (Sfon) ####
# Input trait data, and remove unneeded columns
traitData <- files.df[files.df$species=="Sfon", ]
dim(traitData)
colnames(traitData)
str(traitData)

# Choose phenotypes to retain (choose in order, leave file.name and lib.ID first)
traits.to.keep.all <- c("file.name", "lib.ID" 
                        , "sex", "matur"
                        , "leng.cm_1109", "weight.g_1109", "sp.growth.rateT1.T3"
                        , "condit.fact_T3", "weight_liver.g"
                        , "cort.poststress"
                        , "cort.delta", "osmo.poststress", "osmo.delta"
                        , "chlor.poststress", "chlor.delta"
                        , "fem.egg.diam", "male.sperm.conc", "male.sperm.diam"
                        , "RIN")


allTraits <- traitData[, traits.to.keep.all]
colnames(allTraits)

# Dataframe w/ traits to match expr data
libSamples = rownames(datExpr)
traitRows = match(libSamples, allTraits$file.name) # match w/ file names
datTraits = allTraits[traitRows,] # puts the order of the traits into the order of the samples

rownames(datTraits) <- datTraits[,2] # now use lib.ID as row name
datTraits <- datTraits[,-c(1:2)] # remove file names and lib names
rownames(datTraits)
str(datTraits)
collectGarbage()

# change datTraits to numeric
for(t in 1:length(datTraits)){
  datTraits[,t] <- as.numeric(datTraits[,t])
}
str(datTraits)


# Rename traits
names(datTraits)
names(datTraits) <- gsub(x = names(datTraits), pattern = "RIN", replacement = "RIN")
names(datTraits) <- gsub(x = names(datTraits), pattern = "weight.g_1109", replacement = "weight")
names(datTraits) <- gsub(x = names(datTraits), pattern = "sp.growth.rateT1.T3", replacement = "growth")
names(datTraits) <- gsub(x = names(datTraits), pattern = "leng.cm_1109", replacement = "length")
names(datTraits) <- gsub(x = names(datTraits), pattern = "sex", replacement = "sex")
names(datTraits) <- gsub(x = names(datTraits), pattern = "matur", replacement = "maturity")
names(datTraits) <- gsub(x = names(datTraits), pattern = "condit.fact_T3", replacement = "condition")
names(datTraits) <- gsub(x = names(datTraits), pattern = "weight_liver.g", replacement = "liver weight")
names(datTraits) <- gsub(x = names(datTraits), pattern = "cort.poststress", replacement = "cortisol")
names(datTraits) <- gsub(x = names(datTraits), pattern = "cort.delta", replacement = "cortisol change")
names(datTraits) <- gsub(x = names(datTraits), pattern = "osmo.poststress", replacement = "osmolality")
names(datTraits) <- gsub(x = names(datTraits), pattern = "osmo.delta", replacement = "osmolality change")
names(datTraits) <- gsub(x = names(datTraits), pattern = "chlor.poststress", replacement = "chloride")
names(datTraits) <- gsub(x = names(datTraits), pattern = "chlor.delta", replacement = "chloride change")
names(datTraits) <- gsub(x = names(datTraits), pattern = "fem.egg.diam", replacement = "egg diam.")
names(datTraits) <- gsub(x = names(datTraits), pattern = "male.sperm.conc", replacement = "sperm conc.")
names(datTraits) <- gsub(x = names(datTraits), pattern = "male.sperm.diam", replacement = "sperm diam.")
names(datTraits)


# Rename samples
sample.name.split.list <- strsplit(x = rownames(datExpr), split = "lib", fixed = T )
sample.name.split.mat <- matrix(unlist(sample.name.split.list), ncol=2, byrow=T)
sample.name.rebuild <- gsub(pattern = "_R1", replacement = "", sample.name.split.mat[,2])
# sample.names <- paste("lib", sample.name.rebuild, sep = "")
sample.names <- sample.name.rebuild

# Confirm that they still match
rownames(datExpr)[1:10]
sample.names[1:10]

# and rename the datExpr object
rownames(datExpr) <- sample.names


#### Re-cluster after rem Arctic Charr ####
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = F) # Use color to represent trait values (white = low; red = high; grey = NA)

# Plot sample dendrogram w/ traits
filename <- paste("04_results/", "all_BC", "_samp_clust_and_trait_heatmap.pdf", sep = "")
pdf(file = filename, width = 9, height = 7)

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    cex.dendroLabels = 0.9,
                    main = "",
                    #cex.rowText = 0.3,
                    cex.colorLabels = 0.9,
                    #rowTextAlignment = "left",
                    marAll = c(3,6,1,1)
                    #, abHeight = 180,
                    #abCol = "red"
)

dev.off()