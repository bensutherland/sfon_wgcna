# **WGCNA Analysis for Brook Charr** #

Requirements:
`edgeR` https://bioconductor.org/packages/release/bioc/html/edgeR.html



TODO:    
1. Determine best way to map data for use in WGCNA that can be the same for both species    
2. Test mapping efficiency to either reference genome (substantial difference?)    


Note: this repository is for WGCNA analysis of data prepared using the following repo:
https://github.com/bensutherland/SE-reads_assemble-to-counts.git

At the end of the above repo, read counts will be output from HT-seq, and ready to be normalized and analyzed using WGCNA.

Note that currently the input file has five lines at the end that begin with '__' that are not genes, they are summary information, including 'no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned' and 'alignment_not_unique'. I am not sure what needs to be done with this, as other analysis appears to remove that information.


## 1. Filtering and normalization ##
These steps are performed in edgeR.


