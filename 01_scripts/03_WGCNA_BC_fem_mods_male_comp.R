# wgcna Brook Charr reference module generation and comparison with second dataset

## Clean space
#rm(list=ls())

# Load setup wgcna from previous
load(file = "02_input_data/sfon_wgcna_setup.RData")

#### User: choose working subset ####
## Choose reference set to build modules
# choose among "female", "male" and "AC" for comparisons
REF <- "female"
SEC <- "male"

dim(datExpr.list[[REF]])
datExpr0 <- datExpr.list[[REF]]
dim(datExpr.list[[SEC]]) # this will be input into datExpr0 below during module comparison


#### 2. Filter the subset ####
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
filename <- paste("04_results/", REF, "_transcript_expr_in_num_samp.pdf", sep = "")
pdf(file = filename, width = 5, height = 4)

par(mfrow=c(1,1))
hist(expressed.tally
     , main = "Transcripts expressed in number samples"
     , xlab = "Num. indiv. expr gene"
     , ylab = "Num. transcripts in bin"
     , las = 1) 

dev.off()

# Filter by low expression
datExpr0.filt <- datExpr0[ ,keep.genes]
dim(datExpr0.filt)

# Replace orignal object with filtered
datExpr0 <- datExpr0.filt


#### 3. Data Quality Control ####
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #see tutorial if not true

# plot samples to detect outliers
filename <- paste("04_results/", REF, "_sample_clust_to_detect_outliers.pdf", sep = "")
pdf(file = filename, width = 8.5, height = 4.5)

par(mfrow=c(1,1), mar = c(0,4,2,0), cex = 0.6)
sampleTree <- hclust(dist(datExpr0), method = "average") # default dist metric = euclidean; hclust agglomeration is "average"
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

## Note: Observations: 
## For mature females, 3 main clusters, outlier having large livers

# To improve resolution on trends unrelated to the outlier cluster (i.e. liver weight), cut the tree.
cutline <- 350 # set cutline for females
abline(h = cutline, col = "red") # add to plot

dev.off()


# Identify which samples are in the largest portion of the cut tree
clust = cutreeStatic(sampleTree, cutHeight = cutline, minSize = 10)
table(clust) # clust 1 contains the samples we want to keep
keepSamples = (clust==1)

datExpr = datExpr0[keepSamples, ] # keep only main group
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#### 4.a. Incorporate trait data (Sfon) ####
# Input trait data, and remove unneeded columns
traitData <- files.df[files.df$species=="Sfon", ]
dim(traitData)
colnames(traitData)
str(traitData)

# Choose phenotypes to retain
traits.to.remove <- c(3:5,7,8:10,14,26:27,30,41:45,46) # female
colnames(traitData[, -c(traits.to.remove)]) # confirm traits to keep

# Keep only the desired traits
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

# change datTraits to numeric
for(t in 1:length(datTraits)){
  datTraits[,t] <- as.numeric(datTraits[,t])
}
str(datTraits)


#### 4.b. Re-cluster after rem outliers ####
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = F) # Use color to represent trait values (white = low; red = high; grey = NA)

# Plot sample dendrogram w/ traits
filename <- paste("04_results/", REF, "_samp_clust_and_trait_heatmap.pdf", sep = "")
pdf(file = filename, width = 9, height = 7)

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

dev.off()


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

filename <- paste("02_input_data/", REF, "_sfon_wgcna_save_point_step5.RData", sep = "")
save.image(file = filename)
# load(filename)


#### 5.b. Optional: select only most connected contigs ####
# Characterize connectivity by checking adjacency of all contigs
ADJ = adjacency(datExpr,power=beta1, type="unsigned") #for an unsigned network
#this creates a matrix of gene x genes (Takes a VERY long time, creates large object ADJ)

filename <- paste("02_input_data/", REF, "_sfon_wgcna_save_point_step6.RData", sep = "")
save.image(file = filename)
# load(filename)
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
filename <- paste("02_input_data/", REF, "_sfon_wgcna_save_point_step7.RData", sep = "")
save.image(file = filename)
# load(filename)

#### 5.c. Topological overlap matrix (TOM) ####
# minimize noise and spurious associations by transforming adjacency to Topological Overlap Matrix
# and calculate the corresponding dissimilarity

# Compute the topological overlap matrix based on the adjacency matrix.
dissTOM=TOMdist(ADJ) # default is "unsigned"
save(dissTOM, file ="dissTOM_goes_w_step7")

gc()

# Hierarchical cluster the TOM matrix
hierTOM = hclust(as.dist(dissTOM), method="average")

# Cut tree using dynamic cutting
colorh1 = cutreeDynamic(hierTOM, cutHeight = NULL, minClusterSize = 50)
dynamicColors = labels2colors(colorh1) # Convert numeric lables into colors

# How many contigs were included?
num.transcripts <- length(ADJ[,1])

# Plot dendrogram w/ module colors
filename <- paste("04_results/", REF, "_plotDendroAndColors.pdf", sep = "")
pdf(file = filename, width = 10, height = 5)

par(mfrow=c(1,1))
plotDendroAndColors(hierTOM, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste(c(num.transcripts, "contigs"))
                    , las = 1)
dev.off()

table(dynamicColors) # how many modules were identified and what are the module sizes
unmerged_modules_counts <- as.data.frame(table(dynamicColors)) # how many modules were identified and what are the module sizes
# write.csv(unmerged_modules_counts, file = "04_results/unmerged_modules_counts_fem_filt.csv")

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
filename <- paste("04_results/", REF, "_clust_mod_eigengenes_and_merge_line.pdf", sep = "")
pdf(file = filename, width = 8, height = 5)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=c(0.1,0.2,0.3,0.4), col = c("red", "orange", "green", "blue")) #add level of correlation for cutoff
MEDissThres = 0.25 # note: 0.25 is suggested level from WGCNA Tutorial
abline(h=MEDissThres, col = "purple") # Plot the cut line into the dendrogram

dev.off()

#### 6.b. Merge module eigengenes ####
merge <- mergeCloseModules(datExpr[,restConnectivity], dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors # merged module colors
mergedMEs <- merge$newMEs # merged module eigengenes

# Plot
filename <- paste("04_results/", REF, "_plotDendroAndColors_merged.pdf", sep = "")
pdf(file = filename, width = 10, height = 5)

plotDendroAndColors(hierTOM, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main=paste(num.transcripts,"contigs - merge at",1-MEDissThres))
dev.off()


table(mergedColors) # after merging, how many modules remain and with how many genes
merged_modules_counts <- as.data.frame(table(mergedColors)) # how many modules were identified and what are the module sizes
# write.csv(merged_modules_counts, file = "04_results/merged_modules_counts_0.25_fem_filt.csv")

filename <- paste("02_input_data/", REF, "_sfon_wgcna_save_point_step8.RData", sep = "")
save.image(file = filename)
# load(filename)

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

# Measure dissimilarity b/w module eigengenes (here as a signed correlation)
dissimME <- 1-(t(cor(datMEs, method="pearson")))/2   # spearman is optional if want to try non-parametric

# Cluster based on dissimilarity of module eigengenes
hclustdatME <- hclust(as.dist(dissimME), method="average")

# Plot
filename <- paste("04_results/", REF, "_clust_merged_mod_eigengenes.pdf", sep = "")
pdf(file = filename, width = 8, height = 5)

par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based on the module eigengenes")

dev.off()


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
filename <- paste("04_results/", REF, "_module-trait_rel.pdf", sep = "")
pdf(file = filename, width = 10, height = 9)

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

dev.off()


#### 8 Calculate module membership and gene significance ####
# Per gene/trait combo, Gene Significance *GS* = |cor| b/w  ***gene and trait***
# Per gene/module combo, Module Membership *MM* = cor b/w ***gene and module eigengene***

#### 8.a. Identify trait(s) of interest ####
names(datTraits) # Define variable osmo.delta containing the osmo.delta column of datTrait

# female
TOI.names <- c("weight.g_0709", "sp.growth.rateT1.T3", "condit.fact_T2", "hep.som.ind"
               , "cort.delta", "osmo.delta", "chlor.delta", "fem.egg.diam")

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
#write.csv(geneInfo0, file = "04_results/geneInfo0_fem_filt_most_connected_25000.csv") # FEMALE


#### 9. Identify which female modules are conserved in males ####
#### 9.a. Import male data and filter
datExpr0 <- datExpr.list[[SEC]]

## Filter by defining which genes have > req num expressing samples
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
filename <- paste("04_results/", SEC, "_transcript_expr_in_num_samp.pdf", sep = "")
pdf(file = filename, width = 5, height = 4)

par(mfrow=c(1,1))
hist(expressed.tally
     , main = "Transcripts expressed in number samples"
     , xlab = "Num. indiv. expr gene"
     , ylab = "Num. transcripts in bin"
     , las = 1) 

dev.off()


# Filter by low expression
datExpr0.filt <- datExpr0[ ,keep.genes]
dim(datExpr0.filt)

# Replace orignal object with filtered
datExpr0 <- datExpr0.filt

# Rename the second datExpr0
datExpr0.SEC.low.expr.filt <- datExpr0 

#### 10. Differential network analysis ####

## Samples:
# files.retain.fem
# files.retain.male

#### 10.a. Set up for module preservation ####
# This is the male filtered expression data
dim(datExpr0.SEC.low.expr.filt)
#datExpr.mal <- datExpr0.male.low.expr.filt

## Rename datExpr (from REF, filtered)
datExpr.REF <- datExpr
dim(datExpr.REF)

# Limit to only the top connected genes (previously identified)
datExpr0.REF.mat.top25000 <- datExpr.REF[,restConnectivity]
dim(datExpr0.REF.mat.top25000)

# Set up multi-expression data with module colors from female
setLabels = c(REF, SEC)
multiExpr = list(REF = list(data = datExpr0.REF.mat.top25000), SEC = list(data = datExpr0.SEC.low.expr.filt))
# note: currently SEC is filtered on low expression specifically

goodSamplesGenes(datExpr0.REF.mat.top25000)
goodSamplesGenes(datExpr0.SEC.low.expr.filt)

multiColor = list(REF = mergedColors)

#### 10.b. Test for module preservation ####
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} )


# save(mp, file = "04_results/modulePreservation.RData")
## reload object
# load(file = "04_results/modulePreservation.RData")

#### 10.c. Analyze and visualize module preservation ####
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
# write.csv(mod.pres.table, file = "04_results/mod_perservation_table.csv")


# Set up plot of MedianRank.pres (preservation) & Zsummary.pres for fem mods ~ mod size
plotMods <- !(modColors %in% c("grey", "gold")) # Do not include grey and gold modules
text = modColors[plotMods] # Text labels for points
# For ease of plotting, set up aux variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");

# Plot
filename <- paste("04_results/", REF, "_mod_conservation_in", SEC, ".pdf", sep = "")
pdf(file = filename, width = 10, height = 6)

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

dev.off()


#### 10.d. Plot other statistics in one plot ####
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
filename <- paste("04_results/", REF, "_zstats_in_comp_w_", SEC, ".pdf", sep = "")
pdf(file = filename, width = 12, height = 10)

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

dev.off()

filename <- paste("02_input_data/", REF, "_sfon_wgcna_save_point_step10.RData", sep = "")
save.image(file = filename)
# load(filename)
