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
# # All Brook Charr individuals
# files.retain.BC <- files.df$file.name[files.df$fish.id != "NA"] # no AC
# files.retain.BC
# datExpr0.BC <- datExpr0[files.retain.BC, ]
# dim(datExpr0.BC) # 47 indiv, no parent
# rownames(datExpr0.BC)

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
datExpr0.bck <- datExpr0

# Create list object with all subsets
datExpr.list <- list(female = datExpr0.fem, male = datExpr0.male, AC = datExpr0.AC)
names(datExpr.list)
dim(datExpr.list$female)
dim(datExpr.list[["female"]])

save.image(file = "02_input_data/sfon_wgcna_setup.RData")
