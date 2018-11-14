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
#Make a new variable from maturity and sex
phenodata$sexmatur<-paste(phenodata$sex, phenodata$matur)

# Make a variable indicating outlier samples in network analysis, use in plotting
phenodata$outliers<-phenodata$lib.ID
levels(phenodata$outliers)[levels(phenodata$outliers) %in% c('lib69', 'lib63','lib87', 'lib71','lib85',
                                                          'lib13', 'lib84', 'lib18', 'lib05', 'lib49',
                                                          'lib77', 'lib50', 'lib62')] <- 'yes'
levels(phenodata$outliers)[levels(phenodata$outliers) !='yes'] <- 'no'

# Vector of phenotypes
heterosc<-c("weight_liver.g",
            
            "cort.poststress",
            
            "cort.delta",
            
            "osmo.poststress",
            
            "osmo.delta",
            
            "chlor.delta")

# Vector of y-axis labels
labels<- c("weight_liver.g" = "Liver weight (g)",
           
           "cort.poststress" = expression(paste("Cortisol post-stress (", mu, "g/dL plasma)")),
           
           "cort.delta" = expression(paste("Change in cortisol (", mu ,"g/dL plasma)")),
           
           "osmo.poststress" = "Osmolality post-stress (mmol/kg)",
           
           "osmo.delta" = "Change in osmolality (mmol/kg)",
           
           "chlor.delta" = "Change in chloride (mmol/L)")


## Make box plots for supplemental material

pdf(file="Phenodata_boxplot.pdf", width=7, height=5)
par(mfcol=c(2,3), mar= c(3,5,1,1))
for (x in heterosc) {
  
  # Subset results per trait, write to object 'current.trait'
  current.trait  <- phenodata[,c(x, "outliers", "sexmatur")]
  nonoutliers <- subset(current.trait, outliers == 'no')
  outliers<-subset(current.trait, outliers == 'yes')
  #Boxplot from trait x with sex and maturity
  boxplot(current.trait[,1] ~ current.trait$sexmatur,
          ylab = labels[x])
  #Add points for individuals, first non-outliers
  stripchart(nonoutliers[, 1] ~ nonoutliers$sexmatur, 
            add=T, vertical=T, pch=21, bg="grey", method='jitter', jitter=0.02, cex = 1.1)
  #Add outlier points in blue, indicate locations with 'at' (no outliers in M-).
  stripchart(outliers[,1] ~ outliers$sexmatur, 
             add=T, vertical=T, pch=21, bg="blue", method='jitter', cex = 1.1, at = c(1,2,4))
  
}
  
dev.off()
  