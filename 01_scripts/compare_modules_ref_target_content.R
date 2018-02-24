# Compare modules for content bw sexes

# Input data
# female file: female_geneInfo_25k.csv
# male file: male_geneInfo_25k.csv

# rm(list=ls())

# Set working directory
setwd("~/Documents/10_bernatchez/10_paralogs")


# Import
geneInfo.fem <- read.delim2(file = "female_geneInfo_25k.csv", header = T, sep = ",", row.names = 1)
geneInfo.male <- read.delim2(file = "male_geneInfo_25k.csv", header = T, sep = ",", row.names = 1)


geneInfo.fem  <- geneInfo.fem[,1:3]
geneInfo.male  <- geneInfo.male[,1:3]

head(geneInfo.fem)
head(geneInfo.male)

colnames(geneInfo.fem)
# the colname we are looking for is col3 "moduleColor", and 1 "transcript_id"


fem.mods <- as.character(unique(geneInfo.fem$moduleColor))
male.mods <- as.character(unique(geneInfo.male$moduleColor))

# Which modules are contained within the female modules

##### General Idea pre-looping ####

# Here is the general idea:
head(geneInfo.fem)
mod <- fem.mods[1]

# subset the geneInfo for just one module
subset.fem <- geneInfo.fem[geneInfo.fem$moduleColor==mod, ]
head(subset.fem)
length(subset.fem$transcript_id) # the number of genes in this module

# a variable for the number of genes in this mod
num.tr.in.this.ref.mod <- length(subset.fem$transcript_id)

# How many of the genes in this mod are also in the top 25 k co-expr genes for the second network?
subset.fem.and.male <- merge(x = subset.fem, y = geneInfo.male, by = "transcript_id")
head(subset.fem.and.male)
length(subset.fem.and.male$transcript_id)

# And how many second network modules do these genes belong to in the second network?
length(unique(subset.fem.and.male$moduleColor.y)) # note also can include grey

# Lets summarize how many of this mods genes are in the second network mods
table(subset.fem.and.male$moduleColor.y) # the number of this ref network mod's genes in each non-ref networks mods

# Lets see what percent of the total genes in the ref mod are in each second mod, sorted by descending
sort(round((table(subset.fem.and.male$moduleColor.y) / num.tr.in.this.ref.mod) * 100, digits = 1), decreasing = T)

# Some overviews of total genes in each modules for both networks
sort(table(geneInfo.fem$moduleColor), decreasing = T)
sort(table(geneInfo.male$moduleColor), decreasing = T)

# how many does this represent?
print(mod)

# Lets see for each of the male mods that are represented here, what percentage of that male mod is present in this female mod
sort(table(subset.fem.and.male$moduleColor.y) / table(geneInfo.male$moduleColor) * 100 , decreasing = T)


#### Run in a Loop ####
# First choose sex for reference
#sex <- "fem" ; sec.sex <- "male"
sex <- "male" ; sec.sex <- "fem"


# Then run the following (leave as is)
all.mods.list <- list()
all.mods.list[["fem.mods"]] <- fem.mods
all.mods.list[["male.mods"]] <- male.mods

geneInfo.list <- list()
geneInfo.list[["geneInfo.fem"]] <- geneInfo.fem
geneInfo.list[["geneInfo.male"]] <- geneInfo.male

# The ref mods
ref.mods <- all.mods.list[[paste(sex, "mods", sep = ".")]]
sec.mods <- all.mods.list[[paste(sec.sex, "mods", sep = ".")]]

# The ref geneInfo
ref.geneInfo.name <- paste("geneInfo", sex, sep = ".")
sec.geneInfo.name <- paste("geneInfo", sec.sex, sep = ".")

# set nulls
my.list <- list() ; mod <- NULL ; subset.ref <- NULL ; perc.of.ref <- list() ; subset.ref.and.sec <- NULL; per.missing <- NULL

for(mod in ref.mods){
  print(paste("Next Ref (", sex,") Module Is..."))
  print(mod)
  
  # subset the ref geneInfo data to only the genes w/in the mod of interest
  geneInfo <- geneInfo.list[[ref.geneInfo.name]]
  subset.ref <- geneInfo[geneInfo$moduleColor==mod, ]
  
  # identify the sec geneInfo
  geneInfo.sec <- geneInfo.list[[sec.geneInfo.name]]
  
  # Merge the subset geneInfo with the sec geneInfo for genes found in both
  subset.ref.and.sec <- merge(x = subset.ref, y = geneInfo.sec, by = "transcript_id")
  
  # What percent of the sec modules are represented by the genes found in this ref module?
  my.list[[mod]] <- round(sort(table(subset.ref.and.sec$moduleColor.y) / table(geneInfo.sec$moduleColor) * 100 , decreasing = T), digits = 1)
  # Reporting:
  print(paste("The percentage of each ", sec.sex, "mod that are contained within this", sex, "module", mod, "are:"))
  print(my.list[[mod]])
  
  # What percent of the ref module are found in each of the sec modules here?
  #### ! note this doesn't incl those that aren't in both networks
  num.tr.in.this.ref.mod <- length(subset.ref$transcript_id) 
  perc.of.ref[[mod]] <- sort(round((table(subset.ref.and.sec$moduleColor.y) / num.tr.in.this.ref.mod) * 100, digits = 1), decreasing = T)
  # Reporting:
  print(paste("The percentage of this Ref (", sex, ") module ", mod,  " contained within each ", sec.sex, "modules are:"))
  print(perc.of.ref[[mod]])
  
  # What percent of the ref module is not present in the sec network?
  per.missing <- 100 - sum(perc.of.ref[[mod]])
  print(paste("The % of the ref module that is not present in the non-ref network is", per.missing))
  
}

