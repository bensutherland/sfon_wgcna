## WGCNA analysis
# Second step of the WGCNA repo, with input originating from 01_edgeR_normalization.R    

## Clean space
#rm(list=ls())

## Install Packages
# source("http://bioconductor.org/biocLite.R")
# install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", "fastcluster", "dynamicTreeCut", "survival"))
# biocLite(c("GO.db", "preprocessCore", "impute"))
# install.packages("WGCNA")
require("WGCNA")
#biocLite("edgeR")
require("edgeR")

# Note: currently female data, excluding large liver weight and immature females 
# because manual says to reduce variation and remove outliers

## Set working directory
## Wayne
# setwd("~/Documents/bernatchez/01_Sfon_projects/04_Sfon_eQTL/sfon_wgcna")
# enableWGCNAThreads(nThreads = 3) #Wayne

# Logan & Xavier
setwd("~/Documents/10_bernatchez/01_sfon_eqtl/sfon_wgcna/")
enableWGCNAThreads(nThreads = 7) # Logan
# enableWGCNAThreads(nThreads = 14) # Xavier

# Important setup for loading expression data
options(stringsAsFactors = FALSE) #IMPORTANT SETTING


#### 1.a. Import interp file ####
# Load part 1 results
load("02_input_data/sfon_wgcna_01_output.RData")

### Add Arctic Charr to interp
AC.names <- colnames(normalized.output.log2)[105:122] #contains AC names
# str(interp)

# Generate Arctic char block for interp
number.phenos <- length(colnames(interp))-1 # number phenotypes
number.AC.ind <- length(AC.names) # number of Arctic Charr to add
AC.cells <- number.phenos*number.AC.ind # 792

# Make block dataframe to add
interp.addition <- as.data.frame(matrix(data = c(AC.names, rep("NA", times = AC.cells))
                , nrow=18, ncol=45))
colnames(interp.addition) <- colnames(interp) # make matching column names

# Attach AC to the interp
interp.w.AC <- rbind2(x = interp, y = interp.addition)

### Add species column to the interp and make df
species <- c(rep("Sfon", times = length(interp[,1]))
                 , rep("Salp", times = length(interp.addition[,1])))
interp.final <- cbind(interp.w.AC, species)

files.df <- data.frame(interp.final,stringsAsFactors = F)
# str(files.df)
files.df$sex <- as.character(files.df$sex)
files.df$matur <- as.character(files.df$matur)
str(files.df)
# still all characters... but this is fixed later in "Incorporate trait data.."

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


#### 1.b. Create expression object #### 
sfeqtl <- normalized.output.log2 # (cols = samples, rows = contigs)
dim(sfeqtl)
sfeqtl[1:5,1:5]
str(sfeqtl[1:5,1:5]) #confirm expression values are numeric

# Transpose data
datExpr0 = as.data.frame(t(sfeqtl))
colnames(datExpr0)[1:4] # genes
rownames(datExpr0)[1:4] # samples

# For Testing, Make a Mini Data Set
datExpr0.test <- datExpr0[,1:500]
datExpr0 <- datExpr0.test
dim(datExpr0)

#### 2.a. Create subsets of data (samples) ####
# All Brook Charr individuals
files.retain.BC <- files.df$file.name[files.df$fish.id != "NA"] # no AC
files.retain.BC
datExpr0.BC <- datExpr0[files.retain.BC, ]
dim(datExpr0.BC) # 47 indiv, no parent
rownames(datExpr0.BC)

#females (all, not parent)
files.retain.fem <- files.df$file.name[files.df$sex == "0" & files.df$fish.id != "F2F" 
                                     & files.df$fish.id != "NA"] # Show fish.id, no parent, no AC
files.retain.fem
datExpr0.fem <- datExpr0[files.retain.fem, ]
dim(datExpr0.fem) # 47 indiv, no parent
rownames(datExpr0.fem)

# #males (maturity all same)
files.retain.male <- files.df$file.name[files.df$sex == "1" & files.df$fish.id != "F2M"
                  & files.df$fish.id != "NA"] # Show file name, no parent
files.retain.male
datExpr0.male <- datExpr0[files.retain.male, ]
dim(datExpr0.male) # 53 indiv, no parent
# 

# #Arctic Charr (all males)
salp.interp <- read.csv("02_input_data/salp_interp_table_2017-04-11.csv")
AC.names.8deg <- salp.interp$sample[salp.interp$temp=="8"]

files.retain.AC <- AC.names.8deg
datExpr0.AC <- datExpr0[files.retain.AC, ]
dim(datExpr0.AC) # 10 indiv

# backup all data before subset
#datExpr0.bck <- datExpr0
##### RESTART SUBSET ####
# datExpr0 <- datExpr0.bck

#### 2.b. Choose working subset ####
## Brook Charr all
#datExpr0 <- datExpr0.BC

## female, all
datExpr0 <- datExpr0.fem

## male 
#datExpr0 <- datExpr0.male

## Arctic Charr 8 degrees
#datExpr0 <- datExpr0.AC

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
par(mfrow=c(1,1))
hist(expressed.tally
     , main = "Transcripts expressed in number samples"
     , xlab = "Num. indiv. expr gene"
     , ylab = "Num. transcripts in bin"
     , las = 1) 

# save out as 5 x 4
# transcript_expr_in_num_samples_<sex>.pdf

# Filter by low expression
datExpr0.filt <- datExpr0[ ,keep.genes]
dim(datExpr0.filt)

# Replace orignal object with filtered
datExpr0 <- datExpr0.filt

#### GO ENRICHMENT OUTPUT ####
# # save out as background for GO enrichment
# background <- datExpr0
# annot = read.table(file = "02_input_data/sfontinalis_contigs_annotation_report_v1.0_shortform.txt",
#                    sep = "\t", header = TRUE)
# probes <- names(datExpr0)
# probes2annot <- match(probes, annot$transcript_id) # Index position of probe in annot. file
# sum(is.na(probes2annot)) #number of contigs not present in the annotation file
# tail(probes)
# background <- data.frame(transcript_id = probes,
#                          uniprot_id = annot$sprot_Top_BLASTX_hit[probes2annot])
# head(background)
# 
# # save out background
# write.csv(x = datExpr0, file = "04_results/background_genes_fem_mat.csv")


#### 3. Data Quality Control ####
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #see tutorial if not true

# plot samples to detect outliers
par(mfrow=c(1,1), mar = c(0,4,2,0), cex = 0.6)
sampleTree <- hclust(dist(datExpr0), method = "average") # default dist metric = euclidean; hclust agglomeration is "average"
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# For mature females, 3 main clusters, outlier having large livers
# For males, 2 potential outliers, lib69 and lib37

# To improve resolution on trends unrelated to the outlier cluster (i.e. liver weight), cut the tree.
# Set cutline for males or females
cutline <- 350 # choose height at which to cut (females)
# cutline <- 360 # for males
# cutline <- 600 # to not cut at all (include all samples in trait file)
abline(h = cutline, col = "red") # add to plot

# Save out as 8.5 x 4.5
# sample_clust_to_detect_outliers_<sex>.pdf

clust = cutreeStatic(sampleTree, cutHeight = cutline, minSize = 10)
table(clust) # clust 1 contains the samples we want to keep
keepSamples = (clust==1)

datExpr = datExpr0[keepSamples, ] # keep only main group
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#### FOR MALES ####
# At this point, if working with males save the object to be used as the comparison test set for module conservation
# datExpr0.male.low.expr.filt <- datExpr

#### 4.a. Incorporate trait data (Sfon) ####
# Input trait data, and remove unneeded columns
traitData <- files.df[files.df$species=="Sfon", ]
dim(traitData)
colnames(traitData)
str(traitData)

# Choose phenotypes for females or males
# traits.to.remove <- c(3:5,7,8:10,14,26:27,30,41,43,45,46) # Brook Charr all
traits.to.remove <- c(3:5,7,8:10,14,26:27,30,41:45,46) # female
# traits.to.remove <- c(3:5,7,8:10,14,26:27,30,40:41,43,45,46) # male


colnames(traitData[, -c(traits.to.remove)]) # confirm traits to keep
allTraits <- traitData[, -c(traits.to.remove)] # remove unneeded columns

# Dataframe w/ traits to match expr data
libSamples = rownames(datExpr)
traitRows = match(libSamples, allTraits$file.name) # match w/ file names
datTraits = allTraits[traitRows,] #collect traits for required samples

rownames(datTraits) <- datTraits[,2] # now use lib.ID as row name
datTraits <- datTraits[,-c(1:2)] # remove file names and lib names
rownames(datTraits)
str(datTraits)
collectGarbage()

## backup datTraits object
#datTraits.bck <- datTraits
# datTraits <- datTraits.bck # how to go backwards

# change datTraits to numeric
for(t in 1:length(datTraits)){
  datTraits[,t] <- as.numeric(datTraits[,t])
}
str(datTraits)

# NO LONGER NECESSARY?
# reduce name size for datExpr
# rownames(datExpr) <- sub('.*\\.', '', x = rownames(datExpr)) 
# # this works by searching '.*' one or more characters, followed by a dot '\\.' (\\ is escape), this will match until the last . char in the string.
# rownames(datExpr) <- gsub(x = rownames(datExpr), pattern = "_R1", replacement = "")

# Recap:
# expr data is 'datExpr'
# trait data is 'datTraits'
# and they are matched with lib.id (rownames each)

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
                    marAll = c(3,6,1,1)
                    #, abHeight = 180,
                    #abCol = "red"
                    )

# save as 9 x 7
# sample_clust_and_trait_heatmap_<sex>.pdf

#### 5. Generate Clusters ####
#### 5.a Determine soft-thresh power (Network Topol) ####
# not necessary to re-run after determining best B (soft-thresholding power)
powers = c(1:10, seq(from = 12, to=20, by=2)) # Choose a set of soft-thresholding powers

# ## Call the network topology analysis function
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5
#                        , networkType = "unsigned"  #default
#                        )

# 
# ## Plot to pick soft threshold power
par(mfrow = c(1,2))
cex1 = 0.9
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
# load(file = "02_input_data/sfon_wgcna_save_point_step7_MALE.RData") # load male data


#### 5.c. Topological overlap matrix (TOM) ####
# minimize noise and spurious associations by transforming adjacency to Topological Overlap Matrix
# and calculate the corresponding dissimilarity

# Compute the topological overlap matrix based on the adjacency matrix.
dissTOM=TOMdist(ADJ) # default is "unsigned"
save(dissTOM, file ="dissTOM_goes_w_step7")
# save(dissTOM, file ="dissTOM_goes_w_step7_MALES") # for males
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
# as plotDendroAndColors_<sex>.pdf

table(dynamicColors) # how many modules were identified and what are the module sizes
unmerged_modules_counts <- as.data.frame(table(dynamicColors)) # how many modules were identified and what are the module sizes
write.csv(unmerged_modules_counts, file = "04_results/unmerged_modules_counts_fem_filt.csv")
# write.csv(unmerged_modules_counts, file = "04_results/unmerged_modules_counts_male_filt.csv") #MALES

#### 6. Generate module eigengenes ####
#### 6.a. Create and cluster module eigengenes ####
# Create module eigengenes
MEList <- moduleEigengenes(datExpr[,restConnectivity], colors = dynamicColors)
MEs <- MEList$eigengenes

# Calculate dissimilarity between module eigengenes
MEDiss <- 1-cor(MEs) #dissim

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average") #cluster

# Plot clusters of module eigengenes
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=c(0.1,0.2,0.3,0.4), col = c("red", "orange", "green", "blue")) #add level of correlation for cutoff
MEDissThres = 0.25 # note: 0.25 is suggested level from WGCNA Tutorial
abline(h=MEDissThres, col = "purple") # Plot the cut line into the dendrogram

# save as 8 x 5
# as clust_of_mod_eigengenes_and_merge_line_<sex>.pdf

#### 6.b. Merge module eigengenes ####
merge <- mergeCloseModules(datExpr[,restConnectivity], dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors # merged module colors
mergedMEs <- merge$newMEs # merged module eigengenes

plotDendroAndColors(hierTOM, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main=paste(num.transcripts,"contigs - merge at",1-MEDissThres))
# save out as 10 x 5
# as plotDendroAndColors_merged_<sex>.pdf

table(mergedColors) # after merging, how many modules remain and with how many genes
merged_modules_counts <- as.data.frame(table(mergedColors)) # how many modules were identified and what are the module sizes
write.csv(merged_modules_counts, file = "04_results/merged_modules_counts_0.25_fem_filt.csv")
# write.csv(merged_modules_counts, file = "04_results/merged_modules_counts_0.25_male_filt.csv") #MALE


# save.image(file = "02_input_data/sfon_wgcna_save_point_step8.Rdata")
# save.image(file = "02_input_data/sfon_wgcna_save_point_step8_MALE.Rdata") #MALE

#### 6.c. Correlate module eigengenes
names(mergedMEs)
datMEs <- mergedMEs # datMEs contains module eigengenes' values for each sample

# This assumes the order of the samples is the same as
rownames(datExpr)
### CONFIRM THIS ###

## Export eigengene values e.g. for eQTL of modules
eigengenes.output <- datMEs
rownames(eigengenes.output) <- rownames(datExpr)
dim(eigengenes.output)

# save appropriate dataset
# write.csv(x = eigengenes.output, file = "04_results/eigengenes_output_fem_filt.csv")
# write.csv(x = eigengenes.output, file = "04_results/eigengenes_output_male_filt.csv")


# Measure dissimilarity b/w module eigengenes (here as a signed correlation)
dissimME <- 1-(t(cor(datMEs, method="pearson")))/2   # spearman is optional if want to try non-parametric

# Cluster based on dissimilarity of module eigengenes
hclustdatME <- hclust(as.dist(dissimME), method="average")

# Plot
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based on the module eigengenes")

# save out as 8 x 5
# clust_merged_mod_eigengenes_<sex>.pdf

#### 7.a. Correlate module eigengenes w/ traits and plot ####
# calculate correlation
moduleTraitCor <- cor(mergedMEs, datTraits, use = "p") # p = pairwise.complete.obs (compute in presence of missing values)
# pairwise.complete.obs calcs cor bw ea pair of variables using all complete pairs of obs on those variables

# calculate p-value for correlation
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) #calculates Student asymptotic p-value for given correlations (no MTC)

# Collect results into a matrix for plotting
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #give dimensions of textMatrix

# Plot the correlation values in a heatmap plot
par(mar = c(7, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# fem.mat, 25 modules, 27 traits save out as 10 x 9
# male, 27 modules, 28 traits save out as 10 x 9
# module-trait_rel_<sex>.pdf

#### 8 Calculate module membership and gene significance ####
# Per gene/trait combo, Gene Significance *GS* = |cor| b/w  ***gene and trait***
# Per gene/module combo, Module Membership *MM* = cor b/w ***gene and module eigengene***

#### 8.a. Identify trait(s) of interest ####
names(datTraits) # Define variable osmo.delta containing the osmo.delta column of datTrait

# female
TOI.names <- c("weight.g_0709", "sp.growth.rateT1.T3", "condit.fact_T2", "hep.som.ind"
               , "cort.delta", "osmo.delta", "chlor.delta", "fem.egg.diam")
# male
TOI.names <- c("weight.g_0709", "sp.growth.rateT1.T3", "condit.fact_T2", "hep.som.ind"
               , "cort.delta", "osmo.delta", "chlor.delta", "male.sperm.conc", "male.sperm.diam")



# Create data.frame of traits of interest
TOI = as.data.frame(datTraits[TOI.names])
dim(TOI)

# Obtain names of modules
modNames <- substring(names(mergedMEs), 3) 


#### 8.b. Calculate Module Membership (gene vs eigengene) ####
geneModuleMembership <- as.data.frame(cor(datExpr[,restConnectivity], mergedMEs, use = "p")) #correlate gene with module (MM)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #p-val for MM

# Add names to each dataframe above
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# These are large datasets, and contain values for each transcript against each eigengene
dim(geneModuleMembership)
dim(MMPvalue)


#### 8.c. Calculate Gene Significance (for TOIs) ####
geneTraitSignificance <- as.data.frame(cor(datExpr[,restConnectivity], TOI, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

# Add names to each dataframe above
names(geneTraitSignificance) = paste("GS.", names(TOI), sep="")
names(GSPvalue) = paste("p.GS.", names(TOI), sep="")

# These are large datasets, and contain values for each transcript against each trait
dim(geneTraitSignificance)
dim(GSPvalue)


# #### 8.d. Compare GS and MM for selected module/trait combinations
# 
# names(mergedMEs)
# names(TOI)
# 
# # For example consider
# 
# # note, ME is not present in modNames, and traits have GS. to start
# 
# # female
# module = "firebrick4"
# trait = "GS.cort.delta"
# 
# # male
# module = "steelblue"
# trait = "GS.osmo.delta"
# 
# column <- match(module, modNames) # Index for the module of interest
# column2 <- match(trait, names(geneTraitSignificance)) # Index for trait of interest
# moduleGenes <- mergedColors==module # Index which genes in module of interest
# 
# # Plot MM by GS to see the correlation
# par(mfrow = c(1,1))
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, column2]),
#                    xlab = paste("Module Membership in ", module, " module", sep = ""),
#                    ylab = paste("Gene significance for ", trait, " trait", sep = ""),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

# GS and MM are related
# Theory: Genes highly w/ high GS (gene/trait) oft most important (central) elements of modules assoc with trait


# #### 8.d. Output ####
# # Merge info on modules assoc. w trait of interest (e.g. central genes MM) with gene annotation and write out
# head(names(datExpr)[mergedColors==module]) #just those probes in the module of interest
# length(names(datExpr)[mergedColors==module]) #just those probes in the module of interest

# Read in annotation file
annot = read.table(file = "02_input_data/sfontinalis_contigs_annotation_report_v1.0_shortform.txt",
                   sep = "\t", header = TRUE)
head(annot)
dim(annot)

# Identify probes to annotate
probes = names(datExpr[,restConnectivity])
probes2annot = match(probes, annot$transcript_id) # Index position of probe in annot. file
sum(is.na(probes2annot)) #number of contigs not present in the annotation file
# note that we are currently getting one that is not annotatable, because it is not a real probe "__too_low_aQual" -- should remove these earlier
tail(probes)

# Populate dataframe with key information for all probes
geneInfo0 = data.frame(transcript_id = probes,
                       uniprot_id = annot$sprot_Top_BLASTX_hit[probes2annot],
                       moduleColor = mergedColors,
                       geneTraitSignificance,
                       GSPvalue,
                       geneModuleMembership,
                       MMPvalue)
dim(geneInfo0)

# Write out results
write.csv(geneInfo0, file = "04_results/geneInfo0_fem_filt_most_connected_25000.csv") # FEMALE
# write.csv(geneInfo0, file = "04_results/geneInfo0_male_filt_most_connected_25000.csv") # MALE

#### SAVE OUT COMPLETED MALE DATASET (WITH MODULES GENERATED)
save.image(file = "02_input_data/sfon_wgcna_save_point_step9_MALE.RData")

#### 9. Other sex WGCNA analysis ####
# First to redo all steps using male, return to step (To go back) and then 2.b and proceed through to here again

#### 10. Differential network analysis ####
# Much of this analysis comes from the tutorial found:
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/ModulePreservation/Tutorials/

# data is:
files.retain.fem
files.retain.male




#### 10.a. Calculate module preservation ####
dim(datExpr0.fem) #rows = samples, cols = genes
dim(datExpr0.male.low.expr.filt)

# save out the male low.expr filtered data to re-use later
save(datExpr0.male.low.expr.filt, file = "02_input_data/datExpr0.male.low.expr.filt")


# Old Notes, may not still be conducted..
# In the tutorial, the gene names agree in the two sets, but is not necessary
# Here, we will automatically ignore gene names that cannot be matched across datasets
# but make sure at least half the genes in each module are in common bw ref and test




#### GO BACK AND RE-OBTAIN datExpr as female
load(file = "02_input_data/sfon_wgcna_save_point_step8.Rdata")
# female data
datExpr.fem <- datExpr
dim(datExpr.fem) # should be 32 samples

load(file = "02_input_data/datExpr0.male.low.expr.filt")
datExpr.mal <- datExpr0.male.low.expr.filt
dim(datExpr.mal) # should be 52 samples, and has all genes

# note, since we loaded the step8 save point (from female), the most connected genes will be from female
# need to make sure about this however!!!

# set up the multi-set expression data and corresponding module colors
# As the female data contains the reference modules, carry forward the analysis of the 25000 most connected genes
# and use this as the input
datExpr0.fem.mat.top25000 <- datExpr.fem[,restConnectivity] # restConnectivity gives only top connected, #datExpr0 gives only QC filtered samples
dim(datExpr0.fem.mat.top25000)

setLabels = c("Female", "Male")
multiExpr = list(Female = list(data = datExpr0.fem.mat.top25000), Male = list(data = datExpr0.male.low.expr.filt))
# note: currently male is filtered on low expression specifically

goodSamplesGenes(datExpr0.fem.mat.top25000)
goodSamplesGenes(datExpr0.male.low.expr.filt)

multiColor = list(Female = mergedColors)

# Test for module preservation
?modulePreservation
# For each reference-test pair, the function only uses genes 
# that are in common between the reference and test set. 
# Columns are matched by column names, so column names must be valid.

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} )


save(mp, file = "04_results/modulePreservation.RData")
## reload object
# load(file = "04_results/modulePreservation.RData")

#### 10.b. Analyze and visualize module preservation ####
# Isolate obs stats & Z scores
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

# Isolate mod labels and sizes
modColors <- rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes <- mp$preservation$Z[[ref]][[test]][, 1] # These cap at 1000? Why? #issue# #todo# (also in tutorial)

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

mod.pres.table <- as.data.frame(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
write.csv(mod.pres.table, file = "04_results/mod_perservation_table.csv")


# Set up plot of MedianRank.pres (preservation) & Zsummary.pres for fem mods ~ mod size
plotMods <- !(modColors %in% c("grey", "gold")) # Do not include grey and gold modules
text = modColors[plotMods] # Text labels for points
# For ease of plotting, set up aux variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");

# Plot
par(mfrow = c(1,2), mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust plotting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p],
       col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 7000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4
       , las = 1)
  
  # labeling points does not work with standard tutorial, for first panel
  # option with text labels
  text(x = moduleSizes[plotMods], y = plotData[plotMods, p], labels = text, cex = 0.8, adj = -0.08
  			#, srt = 30 # angle font
  			)
  # option with number labels
  #labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }

# save out as 10 x 6


# Note: There is info in tutorial to inspect whether some modules are driven by outlier samples

#### 10.c. Plot other statistics in one plot ####
# Density and connectivity statistics

# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1]

# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));

# Create numeric labels for each module
nums2match <- seq(1:length(modColors[plotMods]))
cols2match <- modColors[plotMods]
nums.cols.df <- as.data.frame(cbind(nums2match, cols2match))
labs <- nums.cols.df$nums2match

#labs = match(modColors[plotMods], standardColors(50)) # not working, not all cols are in there

# Plot each Z statistic in a separate plot.
par(mfrow = c(4,5), mar = c(3,3,2,1), mgp = c(1.6, 0.4, 0))
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
       xlim = c(20, 2000),
       las = 1)
  text(x = moduleSizes[plotMods], y = statsZ[plotMods, s], labels = labs, cex = 1.3
       , adj = 1.5
       )
  #labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}

# Panel for number and color codes:
plot(x = c(0,10), y = c(0, 10), type = "n"
     , xaxt = 'n'
     , yaxt = 'n'
     , ylab = ""
     , xlab = "Number / Color Key")
cex.legend <- 1.0

# one legend
legend(x = "center", y = "center"
       , legend = paste(nums.cols.df$nums2match, nums.cols.df$cols2match, sep = " ")
       , fill = nums.cols.df$cols2match, cex = cex.legend
       , ncol = 2
       , bty="n"
       , x.intersp=0.4)

# save out as 12 x 10

save.image(file = "02_input_data/completed_step10.RData")


