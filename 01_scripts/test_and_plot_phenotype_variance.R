### Compare variances between sexes across phenotypes
# Jenni Prokkola and Ben Sutherland 2018-11-05

setwd("~/Documents/10_bernatchez/01_sfon_eqtl/sfon_wgcna/")
rm(list=ls())
library(car)

#Read phenotype data
phenodata<-read.csv(file="sfeq_interpretation_v1.csv", header=T)

#Use only the continuous phenotypes included in module comparison available in both sexes
phenodata<-phenodata[, c("lib.ID", "sex", "matur", "weight.g_1109", "leng.cm_1109", "sp.growth.rateT1.T3", 
                         "condit.fact_T3", "weight_liver.g","cort.poststress", "cort.delta", 
                         "osmo.poststress", "osmo.delta", "chlor.poststress", "chlor.delta")]

# Exclude samples (parents) with no phenotypic data
phenodata<-subset(phenodata, weight.g_1109 >0)

# Use Levene's test to compare variance between males and females for each trait.
#Weight and length

boxplot(phenodata$weight.g_1109~phenodata$sex)
leveneTest(weight.g_1109~sex, phenodata)
#Homoscedastic

boxplot(phenodata$leng.cm_1109~phenodata$sex)
leveneTest(leng.cm_1109~sex, phenodata)
#Homoscedastic

# Growth rate T1.T3
boxplot(phenodata$sp.growth.rateT1.T3~phenodata$sex)
leveneTest(sp.growth.rateT1.T3~sex, phenodata)
#Homoscedastic

# Condition factor
boxplot(phenodata$condit.fact_T3~phenodata$sex)
leveneTest(condit.fact_T3~sex, phenodata)
#Homoscedastic

# Liver weight
boxplot(phenodata$weight_liver.g~phenodata$sex)
leveneTest(weight_liver.g~sex, phenodata)
# Heteroscedastic, females have higher variance

# Post-stress cortisol
boxplot(phenodata$cort.poststress~phenodata$sex)
leveneTest(cort.poststress~sex, phenodata)
# Heteroscedastic, females have higher variance

# Change in cortisol
boxplot(phenodata$cort.delta~phenodata$sex)
leveneTest(cort.delta~sex, phenodata)
# Heteroscedastic, females have higher variance

# Osmolality 
boxplot(phenodata$osmo.poststress~phenodata$sex)
leveneTest(osmo.poststress~sex, phenodata)
# Heteroscedastic, females have higher variance

# Change in osmolality 
boxplot(phenodata$osmo.delta~phenodata$sex)
leveneTest(osmo.delta~sex, phenodata)
# Heteroscedastic, females have higher variance

# Chloride 
boxplot(phenodata$chlor.poststress~phenodata$sex)
leveneTest(chlor.poststress~sex, phenodata)
# Homoscedastic

# Change in chloride 
boxplot(phenodata$chlor.delta~phenodata$sex)
leveneTest(chlor.delta~sex, phenodata)
# Heteroscedastic, females have higher variance

# Plot the heteroscedastic phenotypes in relation to maturity status for each sex.
# NOTE only two males were immature.

phenodata$sexmatur<-paste(phenodata$sex, phenodata$matur)

# Vector of phenotypes
heterosc<-c("weight_liver.g",
            
            "cort.poststress",
            
            "cort.delta",
            
            "osmo.poststress",
            
            "osmo.delta",
            
            "chlor.delta")


## Make box plots for supplemental material

par(mfcol=c(2,3), mar= c(2,4,1,1))

for (x in heterosc) {
  
  # Subset results per trait, write to object 'current.trait'
  current.trait  <- phenodata[,c(x, "sexmatur")]
  boxplot(current.trait[,1] ~ current.trait$sexmatur,
          ylab = x)
  stripchart(current.trait[,1] ~ current.trait$sexmatur,
             vertical = TRUE, method = "jitter",
             pch = 21, bg="grey", cex=1.1,
             add = TRUE)
}

