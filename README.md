# **WGCNA Analysis for Brook Charr** #
Note: this repository is for WGCNA analysis of data prepared using the following repo:
https://github.com/bensutherland/SE-reads_assemble-to-counts.git
At the end of the above repo, read counts will be output from HT-seq, and ready to be normalized and analyzed using WGCNA.

Note that currently the input file has five lines at the end that begin with '__' that are not genes, they are summary information, including 'no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned' and 'alignment_not_unique'. I am not sure what needs to be done with this, as other analysis appears to remove that information.

Requirements:
`edgeR` https://bioconductor.org/packages/release/bioc/html/edgeR.html    
`WGCNA` https://cran.r-project.org/web/packages/WGCNA/index.html    

TODO:    
1. Determine best way to map data for use in WGCNA that can be the same for both species    
2. Test mapping efficiency to either reference genome (substantial difference?)    


## 1. Filtering and normalization ##
These steps are performed in edgeR.   
See the script `01_edgeR_normalization.R`   

In brief, transcripts are filtered for low expression (counts > x in y samples), normalized by TMM normalization and library size, and translated to cpm values using the normalized library size.
Data is visualized using mds plot, dispersions estimated (coefficient of variation for experiment).
Data is output in a matrix of filtered and normalized genes (but all are still in linear, not log2).
Output datafile can be found in `03_normalized_data/normalized_output_matrix.csv`

Note: original analysis separately generated male and female output and filtering. I don't think this is necessary, as in the next steps, the data can be subset by sex and re-filtered for low expression. The only place that may be negatively impacted is by using the TMM normalization on the sexes together, as we can expect there to be expression differences between these groups. For now will leave as is.

## 2. WGCNA module building in females ##
These steps are performed in R using the WGCNA package.


