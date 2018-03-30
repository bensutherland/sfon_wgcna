#########################################################################
### 
### R code for heatmaps of modules from WGCNA (Sutherland et al. manuscript)
### Jenni Prokkola 
###
#########################################################################
library(RColorBrewer)
library(gplots)
library(doBy)
library(data.table)

#######################################################################################
### Import all required files: 
#######################################################################################
## Female modules
Female_mod<-read.csv(file="female_GeneInfo_25k.csv", header=T, dec=".")
row.names(Female_mod)<-Female_mod[,1]
Female_mod<-Female_mod[,-1]
str(Female_mod)

### Male modules
Male_mod<-read.csv(file="male_GeneInfo_25k.csv", header=T, dec=".")
row.names(Male_mod)<-Male_mod[,1]
Male_mod<-Male_mod[,-1]
str(Male_mod)

### Sex and library information
Sampledata<-read.csv(file="sfeq_interpretation_v1.csv", header=T)
Lib_IDs_reads<-read.csv2(file="Sfon_lib_IDs.csv", header=T)

### Counts data
Log_cpm_data<-read.csv(file="normalized_output_matrix_log2.csv", header=T, dec=".")
### Make gene names row names
row.names(Log_cpm_data)<-Log_cpm_data[,1]
Log_cpm_data<-Log_cpm_data[,-1]

### Edit counts data to use short library IDs
Log_cpm_data_sort<-Log_cpm_data
names(Log_cpm_data_sort)[1:104]<-unlist(strsplit(names(Log_cpm_data_sort)[1:104], "_"), use.names = FALSE)[c(FALSE, TRUE, FALSE)]
names(Log_cpm_data_sort)[1:104]<-unlist(strsplit(names(Log_cpm_data_sort)[1:104], ".", fixed=TRUE), use.names = FALSE)[c(FALSE, TRUE)]
str(Log_cpm_data_sort)

### Combine library ID and sex info
Sex_info<-data.frame(Lib_id=Sampledata$lib.ID, sex= Sampledata$sex)
Lib_IDs_reads<-merge(Lib_IDs_reads, Sex_info, by.x="Lib_ID", by.y="Lib_id", sort=F)
## Add Arctic char males (M2)
Lib_IDs_reads<-rbind(Lib_IDs_reads, data.frame(Lib_ID = names(Log_cpm_data)[105:113], sex = rep("M2", 9)))

## Order the file by sex
Lib_IDs_ordered<-orderBy(~sex, data=Lib_IDs_reads)

## Order counts data based on sex
setDT(Log_cpm_data_sort)
setcolorder(Log_cpm_data_sort, as.character(Lib_IDs_ordered$Lib_ID))
Log_cpm_data_sort<-as.data.frame(Log_cpm_data_sort)

## add gene names again
rownames(Log_cpm_data_sort)<-row.names(Log_cpm_data)

## Make a vector of colors based on sex
x<-as.character(Lib_IDs_ordered$sex)
x[which(x=="M")]<-"steelblue4"
x[which(x=="F")]<-"orange"
x[which(x=="M2")]<-"yellow"

#######################################################################################
## Generating heatmap from single modules, both columns and rows clustered by Pearson correlation distances
## Colors from RColorBrewer
#######################################################################################
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

### Pick module to visualize (male or female data, change module color)
selected_clust<-row.names(subset(Male_mod, moduleColor=="yellow"))
## how many genes?
length(selected_clust)

### Subset the selected genes from counts data
heat_data<- as.matrix(Log_cpm_data_sort[selected_clust,])

### Distances of genes and samples by Pearson correlation
hc_heat_data_row <- hclust(as.dist(1-cor(t(heat_data))), method="average") 
hc_heat_data_col <- hclust(as.dist(1-cor(heat_data)), method="average")

### Draw heatmap to pdf
pdf(file="heatmap_male_yellow.pdf")
heatmap.2(heat_data,Rowv = as.dendrogram(hc_heat_data_row), scale="row", 
          Colv = as.dendrogram(hc_heat_data_col),
          col=hmcol, trace="none", labRow=NA,labCol = NA,
          ColSideColors=c(x),
          margins= c(8,5), srtCol = 70,
          cexCol=0.9, key=TRUE, density.info="density", symkey=FALSE, 
          key.par=list(cex.axis=0.8), key.title=NA, keysize=1)

dev.off()


#######################################################################################################
## Generating heatmap from two modules, both columns and rows clustered by Pearson correlation distances
#######################################################################################################
selected_clust_2<-row.names(rbind(subset(Male_mod, moduleColor=="darkmagenta"), subset(Male_mod, moduleColor=="steelblue")))

### Subset the selected genes from counts data
heat_data<- as.matrix(Log_cpm_data_sort[selected_clust_2,])

### Distances of genes and samples by Pearson correlation
hc_heat_data_row <- hclust(as.dist(1-cor(t(heat_data))), method="average") 
hc_heat_data_col <- hclust(as.dist(1-cor(heat_data)), method="average")

### Draw heatmap to pdf
pdf(file="heatmap_male_darkmagenta_steelblue.pdf")
heatmap.2(heat_data,Rowv = as.dendrogram(hc_heat_data_row), scale="row", 
          Colv = as.dendrogram(hc_heat_data_col),
          RowSideColors=c(rep("darkmagenta", 61), rep("steelblue", 65)),
          col=hmcol, trace="none", labRow=NA,labCol = NA,
          ColSideColors=c(x),
          margins= c(8,5), srtCol = 70,
          cexCol=0.9, key=TRUE, density.info="density", symkey=FALSE, 
          key.par=list(cex.axis=0.8), key.title=NA, keysize=1)

dev.off()



