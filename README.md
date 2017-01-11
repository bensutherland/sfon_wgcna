# **WGCNA Analysis for Brook Charr** #

Requirements:
`edgeR` https://bioconductor.org/packages/release/bioc/html/edgeR.html



TODO:    
1. Determine best way to map data for use in WGCNA that can be the same for both species    
2. Test mapping efficiency to either reference genome (substantial difference?)    


Note: this repository is for WGCNA analysis of data prepared using the following repo:
https://github.com/bensutherland/SE-reads_assemble-to-counts.git

At the end of the above repo, read counts will be output from HT-seq, and ready to be normalized and analyzed using WGCNA.


## 1. Normalization ##
Normalization is performed using edgeR
