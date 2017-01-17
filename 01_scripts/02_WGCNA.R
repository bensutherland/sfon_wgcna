## WGCNA analysis
# This is the second step of the WGCNA repo, with input originating from 01_edgeR_normalization.R    

# rm(list=ls())

## Install Packages
# source("http://bioconductor.org/biocLite.R")
# install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", "fastcluster", "dynamicTreeCut", "survival"))
# biocLite(c("GO.db", "preprocessCore", "impute"))
# install.packages("WGCNA")
require("WGCNA")
#biocLite("edgeR")
require("edgeR")

# Note, the current analysis begins with female data that excludes immature females and large liver weight females, b/c manual says reduce variation and rem. outliers

# Set working directory
setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/sfon_wgcna")
# setwd("~/Documents/sfon_wgcna/") # macpro

#### 1 Import interp file and data ####
# Important setup for loading expression data
options(stringsAsFactors = FALSE) #IMPORTANT SETTING

# Load part 1 results
load("02_input_data/sfon_wgcna_01_output.RData")

# Enable parallel processing
enableWGCNAThreads(nThreads = 2)
# enableWGCNAThreads(nThreads = 10) #MacPro

### End front matter ###

# Create data.frame
files.df <- as.data.frame(interp)
names(files.df)

# There were some issues with the sex of individuals in the interp file, but it looks like those samples are not in this set
files.df$fish.id[files.df$male.sperm.conc != "NA" & files.df$sex == "F"]
files.df$fish.id[files.df$fem.egg.diam != "NA" & files.df$sex == "M"]

# Recode sex as binary
head(files.df[,c("sex","matur")])
files.df$sex[files.df$sex=="F"] <- 0
files.df$sex[files.df$sex=="M"] <- 1
files.df$sex <- as.numeric(files.df$sex)

# Recode matur as binary
files.df$matur[files.df$matur=="-"] <- 0
files.df$matur[files.df$matur=="+"] <- 1
files.df$matur <- as.numeric(files.df$matur)
head(files.df[,c("sex","matur")])

# Make object with expression data
sfeqtl <- normalized.output.log2 # (cols = samples, rows = contigs)
dim(sfeqtl)
sfeqtl[1:5,1:5]
str(sfeqtl[1:5,1:5]) #confirm expression values are numeric

# Transpose data
datExpr0 = as.data.frame(t(sfeqtl))
colnames(datExpr0)[1:4] # genes
rownames(datExpr0)[1:4] # samples

datExpr0.bck <- datExpr0
## To go back
# datExpr0 <- datExpr0.bck


#### 2.a. Create subsets of data (samples) ####
# Change files.df$file.name to match the short form name
files.df$file.name <- gsub(x = files.df$file.name, pattern = ".txt", replacement = "")
head(files.df$file.name)

#females (all, not parent)
files.df$fish.id[files.df$sex == "0" & files.df$fish.id != "F2F"] # Show fish.id, no parent
files.retain.fem <- files.df$file.name[files.df$sex == "0" & files.df$fish.id != "F2F"] # get filenames for subset
files.retain.fem
datExpr0.fem <- datExpr0[files.retain.fem, ]
dim(datExpr0.fem) # 47 indiv, no parent
rownames(datExpr0.fem)

#females (only mature)
files.df$fish.id[files.df$sex == "0" & files.df$matur == 1 & files.df$fish.id != "F2F"] # Show fish.id, no parent
files.retain.fem.mat <- files.df$file.name[files.df$sex == "0" & files.df$matur == 1 & files.df$fish.id != "F2F"]  # get filenames for subset
datExpr0.fem.mat <- datExpr0[files.retain.fem.mat,]
dim(datExpr0.fem.mat) # 41 indiv, no parents
rownames(datExpr0.fem.mat)

#males (maturity all same)
files.df$fish.id[files.df$sex == "1" & files.df$fish.id != "F2M"] # Show fish.id, no parent
files.retain.male <- files.df$file.name[files.df$sex == "1" & files.df$fish.id != "F2M"] # get filenames for subset
datExpr0.male <- datExpr0[files.retain.male, ]
dim(datExpr0.male) # 53 indiv, no parent

#### 2.b. Choose working subset and filter ####
## female, mature
datExpr0 <- datExpr0.fem.mat

# Note: Will need to do some filtering similar to the following 
# (remember, at log2), and cpm thresh was 0.5, so log2(0.5)

#### 2.c. Filter the subset ####
## Define which genes have > req num expressing samples
num.indiv <- 5

dim(datExpr0)
expressed <- NULL
keep.genes <- NULL
expressed.tally <- NULL
for(i in 1:length(colnames(datExpr0))) { 
  expressed <- length(which(datExpr0[,i] > -1))
  print(expressed)
  expressed.tally <- c(expressed.tally, expressed)
  if (expressed > num.indiv) {
    keep.genes <- c(keep.genes, colnames(datExpr0)[i])
  }
}

head(keep.genes)
length(keep.genes)

# Plot the freq of genes for each bin of number samples expressed?
hist(expressed.tally
     , main = "Transcripts expressed in number samples"
     , xlab = "Num. indiv. expr gene"
     , ylab = "Num. transcripts in bin"
     , las = 1) 

# Filter by low expression
datExpr0.filt <- datExpr0[ ,keep.genes]
dim(datExpr0.filt)

# Replace orignal object with filtered
datExpr0 <- datExpr0.filt


#### 3. Data Quality Control ####
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #see tutorial if not true

# plot samples to detect outliers
par(mfrow=c(1,1), mar = c(0,4,2,0), cex = 0.6)
sampleTree <- hclust(dist(datExpr0), method = "average") # default dist metric = euclidean; hclust agglomeration is "average"
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# For mature females, 3 main clusters, outlier having large livers
# To improve resolution on trends unrelated to the outlier cluster (i.e. liver weight), cut the tree.

# Remove Outliers
cutline <- 200 # choose height at which to cut
abline(h = cutline, col = "red") # add to plot

clust = cutreeStatic(sampleTree, cutHeight = cutline, minSize = 10)
table(clust) # clust 1 contains the samples we want to keep
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ] # keep only main group
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#### 4.a. Incorporate trait data ####
# Input trait data, and remove unneeded columns
traitData <- files.df
dim(traitData)
colnames(traitData)
traits.to.remove <- c(3:5,7,8:10,14,26:27,30,41:45)
colnames(traitData[, -c(traits.to.remove)]) # these are the phenos to remove
allTraits <- traitData[, -c(traits.to.remove)] # remove unneeded columns

# Dataframe w/ clinical traits to match expr data
libSamples = rownames(datExpr)
traitRows = match(libSamples, allTraits$file.name)
datTraits = allTraits[traitRows,] #collect traits for required samples
rownames(datTraits) <- datTraits[,2] # use lib.ID as row name
datTraits <- datTraits[,-c(1:2)] # remove file names and lib names
rownames(datTraits)
collectGarbage()


# Recap:
# expr data is 'datExpr'
# trait data is 'datTraits'

# reduce name size for datExpr
rownames(datExpr) <- sub('.*\\.', '', x = rownames(datExpr)) 
# this works by searching '.*' one or more characters, followed by a dot '\\.' (\\ is escape), this will match until the last . char in the string.

#### 4.b. Re-cluster after rem outliers ####
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = F) # Use color to represent trait values (white = low; red = high; grey = NA)

# Plot sample dendrogram w/ traits
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    cex.dendroLabels = 0.7,
                    main = "",
                    #cex.rowText = 0.3,
                    cex.colorLabels = 0.7,
                    #rowTextAlignment = "left",
                    marAll = c(3,6,1,1),
                    #abHeight = 180,
                    #abCol = "red"
                    )


#### 5.a Determine soft-thresh power (Network Topol) ####
# not necessary to re-run after determining best B (soft-thresholding power)
powers = c(1:10, seq(from = 12, to=20, by=2)) # Choose a set of soft-thresholding powers

# ## Call the network topology analysis function
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5
#                         , networkType = "unsigned"  #default 
#                         , ) 
# 
# ## Plot to pick soft threshold power
# par(mfrow = c(1,2))
# cex1 = 0.9
# 
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# 
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Result summary:
# Females, no outliers, choose power = 6 (note that this is the default)
# this was chosen because...

beta1=6 # If using signed network, double the beta1


#### SAVE POINT ####
# save.image(file = "02_input_data/sfon_wgcna_save_point_step5.RData")
# load("02_input_data/sfon_wgcna_save_point_step5.RData")

#### 5.b. Optional: select only most connected contigs ####
# Characterize connectivity by checking adjacency of all contigs
ADJ = adjacency(datExpr,power=beta1, type="unsigned") #for an unsigned network
#this creates a matrix of gene x genes (Takes a VERY long time, creates large object ADJ)

# save.image(file = "02_input_data/sfon_wgcna_save_point_step6.RData")
# load("02_input_data/sfon_wgcna_save_point_step6.RData")
gc()

#calc the connectivity of each gene
Connectivity=softConnectivity(datExpr,power=beta1)-1
gc()

# Restrict to most connected genes
ConnectivityCut = 25000 # number of most connected genes to keep
ConnectivityRank = rank(-Connectivity)
restConnectivity = ConnectivityRank <= ConnectivityCut # true/false for ea. gene to keep or not
sum(restConnectivity) # should match the ConnectivityCut above
gc()

# Re-define adjacency matrix for only those most connected genes
ADJ = adjacency(datExpr[,restConnectivity], power=beta1, type = "unsigned")
gc()

## This should be smaller than the above (4 Gb vs 12 Gb), so re-save out
# save.image(file = "02_input_data/sfon_wgcna_save_point_step7.RData")


#### 5.c. Topological overlap matrix (TOM) ####
# minimize noise and spurious associations by transforming adjacency to Topological Overlap Matrix
# and calculate the corresponding dissimilarity

# Compute the topological overlap matrix based on the adjacency matrix.
dissTOM=TOMdist(ADJ) # default is "unsigned"
gc()

# Hierarchical cluster the TOM matrix
hierTOM = hclust(as.dist(dissTOM), method="average")

# Cut tree using dynamic cutting
colorh1 = cutreeDynamic(hierTOM, cutHeight = NULL, minClusterSize = 50)
dynamicColors = labels2colors(colorh1) # Convert numeric lables into colors

# How many contigs were included?
num.transcripts <- length(ADJ[,1])

# Plot dendrogram w/ module colors
par(mfrow=c(1,1))
plotDendroAndColors(hierTOM, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste(c(num.transcripts, "contigs"))
                    , las = 1)
# save out as 10 x 5

table(dynamicColors) # how many modules were identified and what are the module sizes

#### 5.d. Module eigengenes ####
# Create module eigengenes
MEList <- moduleEigengenes(datExpr[,restConnectivity], colors = dynamicColors) # may be an issue here, bc dynamic colors is on top 25000 genes
MEs <- MEList$eigengenes

# Calc module eigengene dissimilarity matrix
MEDiss <- 1-cor(MEs) #dissim
METree <- hclust(as.dist(MEDiss), method = "average") #cluster

# Plot module eigengene clustering
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=c(0.1,0.2,0.3,0.4), col = c("green", "pink", "blue", "red")) #add level of correlation for cutoff
MEDissThres = 0.25 # 0.25 is suggested level from WGCNA Tutorial
abline(h=MEDissThres, col = "green") # Plot the cut line into the dendrogram

# Merge eigengenes
# top genes:
merge <- mergeCloseModules(datExpr[,restConnectivity], dynamicColors, cutHeight = MEDissThres, verbose = 3)
# merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) # SAME ISSUE AS ABOVE
mergedColors <- merge$colors # merged module colors
mergedMEs <- merge$newMEs # merged module eigengenes

plotDendroAndColors(hierTOM, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main=paste(num.transcripts,"contigs - merge at",1-MEDissThres))
# save out as 10 x 5

table(mergedColors) # after merging, how many modules remain and with how many genes

save.image(file = "02_input_data/sfon_wgcna_save_point_step8.Rdata")



#### FRONTLINE #####

########2C CORRELATE MODULE EIGENGENES WITH EACH OTHER####
names(mergedMEs)
datMEs <- mergedMEs #rename to use in the following

#datMEs contains module eigengenes for each sample
# currently assuming that the order of the samples is the same as:
rownames(datExpr)

#####TEMP ####
##### EQTL OUTPUT ######
eigengenes.output <- datMEs
rownames(eigengenes.output) <- rownames(datExpr)
dim(eigengenes.output)
write.csv(x = eigengenes.output, file = "eigengenes.output.10cpm.femclean.csv")
?write.csv




# Dissimilarity measure bw module eigengenes (signed correlation)
dissimME <- 1-(t(cor(datMEs, method="pearson")))/2   #could use spearman if want to try non-parametric
hclustdatME <- hclust(as.dist(dissimME), method="average")
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based on the module eigengenes")

########3A CORRELATE MODULE EIGENGENES WITH TRAITS ####
moduleTraitCor <- cor(mergedMEs, datTraits, use = "p"); #use=p =pairwise.complete.obs; method for computing covariances in the presence of missing values
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples); #calculates Student asymptotic p-value for given correlations (no MTC)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #give dimensions of textMatrix
par(mar = c(4, 7, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

########3B CORRELATE GENE TO TRAIT AND KEY MODULES (gene sig & mod memb) ####
# Per gene, define Gene Significance *GS* as: the |corr| bw the ***gene and the trait***
# Per gene[module], define Module Membership *MM* as: the corr bw the ***gx profile and the module eigengene***

# identify and separate trait(s) of interest
names(datTraits) # Define variable osmo.delta containing the osmo.delta column of datTrait
TOI = as.data.frame(datTraits[c(8, 19, 22, 25)])

modNames <- substring(names(MEs), 3) # obtain names (colors) of the modules

#Module Membership
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p")) #correlate gene with module (MM)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #p-val for MM
names(geneModuleMembership) = paste("MM", modNames, sep="") #adds names to dataframe
names(MMPvalue) = paste("p.MM", modNames, sep="") #adds names to dataframe

#Gene Significance (for interesting traits)
geneTraitSignificance <- as.data.frame(cor(datExpr, TOI, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(TOI), sep="")
names(GSPvalue) = paste("p.GS.", names(TOI), sep="")

########3C INTRAMODULAR ANALYSIS: GENES W HIGH GS & MM ####
# As an example, we look at the 'white' module that has the highest association with 'osmo.delta' 
module = "white"
trait = "GS.osmo.delta"
column <- match(module, modNames) #column index for the module of interest
column2 <- match(trait, names(geneTraitSignificance)) #column index for trait of interest
moduleGenes <- mergedColors==module #indexes genes that are in the module of interest

# Plot MM by GS to see the correlation
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, column2]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for osmo.delta",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
#Clearly, GS and MM are highly correlated
# Genes highly signif. assoc. w a trait oft also most 
#  important (central) elements of modules assoc with trait

########3D SUMMARY AND OUTPUT (ANNOTATION OF RESULTS) ####
#Merge info on modules assoc. w trait of interest (e.g. central genes MM)
# with gene annotation and write out

# read in annotation file
head(names(datExpr)[mergedColors=="white"]) #just those probes in the module of interest

annot = read.table(file = "/Users/wayne/Documents/bernatchez/Sfon_projects/Sfon-transcriptome/annotating_transcriptome/sfontinalis_contigs_annotation_report_v1.0_shortform.txt",
                   sep = "\t", header = TRUE)
head(annot)
dim(annot)

# identify probes to annotate
probes = names(datExpr)
probes2annot = match(probes, annot$transcript_id) #informs on the positions of the probe in the annotation file
sum(is.na(probes2annot)) #number of contigs not present in the annotation file

# populate dataframe with key information for all probes
geneInfo0 = data.frame(transcript_id = probes,
                       uniprot_id = annot$sprot_Top_BLASTX_hit[probes2annot],
                       moduleColor = mergedColors,
                       geneTraitSignificance,
                       GSPvalue)

osmo.delta <- datTraits$osmo.delta

# order modules by their significance for osmo.delta
modOrder = order(-abs(cor(MEs, osmo.delta, use = "p")))

# add in module membership information in the above order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

dim(geneInfo0)

# write out results
write.csv(geneInfo0, file = "geneInfo0.csv")

########4 DIFFERENTIAL NETWORK ANALYSIS########
########4A MALE EXPRESSION DATA INPUT #######
# use the same files.df file as above, but select only males
mal.files <- files.df$file.name[files.df$sex==1 & files.df$matur==1]
mal.files <- mal.files[!is.na(mal.files)] # this removes parents
posit.mal.files <- which(files.df$file.name %in% mal.files)
mal.files.df <- files.df[posit.mal.files,]
dim(mal.files.df)

# Read in the normalized data (cols = samples, rows = contigs)
# Input male-specific filtered and normalized genes:
sfeqtl <- read.csv(file="/Users/wayne/Documents/bernatchez/Sfon_projects/SfeQ/02-RNA-Seq_analysis/07_gx_levels-may25-15/log2.outcpm.mal_only-cpm5_in25.csv")
dim(sfeqtl)
head(sfeqtl[1:5,1:5])
str(sfeqtl[1:5,1:5]) #confirm variables are numeric

datExpr0 = as.data.frame(t(sfeqtl[, -c(1)])) #remove auxillary data and transpose (excluding contigs)
names(datExpr0) = sfeqtl[,1] #provide contig IDs as colnames
head(datExpr0[1:5,1:5])

# Subset to keep the required male expression data
datExpr0.mal <- datExpr0[mal.files,]
dim(datExpr0.mal)
head(datExpr0.mal[1:5,1:5])

# Replace the long form row.names with short form (libID)
posit.mal.files <- which(mal.files.df$file.name %in% mal.files)
rownames(datExpr0.mal) <- mal.files.df$lib.ID[posit.mal.files]
head(datExpr0.mal[1:5,1:5])

#########4B MALE Expression data checking#####
gsg = goodSamplesGenes(datExpr0.mal, verbose = 3)
gsg$allOK #see tutorial if not true

# plot samples to detect outliers
par(mfrow=c(1,1))
sampleTree.mal <- hclust(dist(datExpr0.mal), method = "average") # default dist metric = euclidean; hclust agglomeration is "average"
par(mar = c(0,4,2,0), cex = 0.6)
plot(sampleTree.mal, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#definitely an odd sample 69, cut at 100

## To improve resolution on trends unrelated to the outlier cluster (i.e. liver weight), cut the tree.
cutline <- 100 # choose height at which to cut
abline(h = cutline, col = "red") # add to plot
clust = cutreeStatic(sampleTree.mal, cutHeight = cutline, minSize = 10)
table(clust) # clust 1 contains the samples we want to keep
keepSamples = (clust==1)
datExpr.mal = datExpr0.mal[keepSamples, ] # subset the expression data to only keep samples in the main group
nGenes = ncol(datExpr.mal)
nSamples = nrow(datExpr.mal)

# expression without outliers is now in datExpr.mal

##########4C MALE TRAIT DATA#######
# can use the original 'allTraits' object as this contais all samples
names(allTraits)

# Form a dataframe analogous to expression data that will hold the clinical traits.
libSamples = rownames(datExpr.mal)
traitRows = match(libSamples, allTraits$lib.ID); #matches sample names! perfect
datTraits.mal = allTraits[traitRows, -1] #collect traits for required samples, remove sample names
rownames(datTraits.mal) = allTraits[traitRows, 1]; #add sample names as row-names
rownames(datTraits.mal)
collectGarbage()

# Post-outlier removal, Re-cluster samples
sampleTree2 = hclust(dist(datExpr.mal), method = "average") 
datTraits.mal.trimmed <- datTraits.mal[,-29] # has trouble with the NAs in the egg category 
#ISSUE HERE ##### NEED TO RETAIN THE SPERM VALUES####
traitColors = numbers2colors(datTraits.mal.trimmed, signed = F) # Use color to represent trait values (white = low; red = high; grey = NA)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits.mal),
                    main = "Sample dendrogram and trait heatmap")

##  expr data is 'datExpr.mal'
#   trait data is 'datTraits.mal'

#########4D CALCULATION OF MODULE PRESERVATION#######
dim(datExpr.mal) #rows = samples, cols = genes
dim(datTraits.mal)

# in the tutorial, the gene names agree in the two sets, but is not necessary
# here, we will automatically ignore gene names that cannot be matched across datasets
# but make sure at least half the genes in each module are in common bw ref and test

#set up the multi-set expression data and corresponding module colors
setLabels = c("Female", "Male")
multiExpr = list(Female = list(data = datExpr), Male = list(data = datExpr.mal))
multiColor = list(Female = mergedColors)
?modulePreservation
#For each reference-test pair, the function only uses genes (columns of the data component of each component of multiExpr) that are in common between the reference and test set. Columns are matched by column names, so column names must be valid.
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} )


save(mp, file = "modulePreservation.RData")

########4E Analyze and Display module Preservation Results#####
#load(file = "modulePreservation.RData")
# isolate the observed statistics and their Z scores
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# plot the preservation medianRank and Zsummary for the female modules as a function of module size

# Module labels and module sizes are also contained in the results
modColors <- rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes <- mp$preservation$Z[[ref]][[test]][, 1] # Why are these capping at 1000? is it just coincidence...??
# leave grey and gold modules out
plotMods <- !(modColors %in% c("grey", "gold"))
# Text labels for points
text = modColors[plotMods]
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");

par(mfrow = c(1,2), mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }

#THERE IS SOMETHING WRONG HERE BECAUSE MY TEXT AREN'T SHOWING UP

#without loop
par(mfrow = c(1,2), mar = c(4.5,4.5,2.5,1))
plot(moduleSizes[plotMods], plotData[plotMods, 1], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[1],
       cex = 2.4,
       ylab = mains[1], xlab = "Module size", log = "x",
       ylim = c(0,15),
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(moduleSizes[plotMods], plotData[plotMods, 1], text, cex = 1, offs = 0.08);

plot(moduleSizes[plotMods], plotData[plotMods, 2], col = 2, bg = modColors[plotMods], pch = 21,
     main = mains[2],
     cex = 2.4,
     ylab = mains[2], xlab = "Module size", log = "x",
     ylim = c(0,40),
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(moduleSizes[plotMods], plotData[plotMods, 2], text, cex = 1, offs = 0.08);
# For Zsummary, add threshold lines
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)


#Closer inspection of the modules with low preservation may indicate that the 
# module was driven by outlier sample (see Supplementary Text S5) in Langfelder et al 2011

#now plot the density and connectivity statistics all in one plot

# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
par(mfrow = c(4,5), mar = c(3,3,2,1), mgp = c(1.6, 0.4, 0))
# Plot each Z statistic in a separate plot.
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/5) min = -max/5
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.7,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(20, 1000))
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}

data.frame(color = modColors[plotMods], label = labs)






