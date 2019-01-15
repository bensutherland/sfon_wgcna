# **WGCNA Analysis for Brook Charr** #
WGCNA analysis of transcriptome data for the manuscript Sex-specific co-expression networks and sex-biased gene expression in the salmonid Brook Charr Salvelinus fontinalis (https://doi.org/10.1101/305680).   
This pipeline has only been tested and used for the purposes of the authors of the above article, and has no guarantees of usefulness whatsoever for other purposes.   
Section B below, which focus on module-chromosome analysis, was used briefly in the above article for analysis of the sex chromosome only.     

*Required inputs to the repo include:*   
1) a counts file with linear values that are output from the following repo:  
https://github.com/bensutherland/Simple_reads_to_counts    
specifically, within this other repo using `01_script/03_express.sh` and then `01_scripts/utility_scripts/prepare_gxlevels_matrix.R`.   
Essentially, using the 'effective counts' column as suggested by eXpress (Roberts and Pachter 2013).     
Put this output into the present repo as: `02_input_data/out.matrix.csv`    

2) an interpretation file that contains the filenames for each sample and the associated phenotypic data. 
This should be entitled/located as follows: `00_archive/sfeq_interpretation_v1.csv`   

*Required libraries:*  
`edgeR` https://bioconductor.org/packages/release/bioc/html/edgeR.html    
`WGCNA` https://cran.r-project.org/web/packages/WGCNA/index.html    
`locfit` https://cran.r-project.org/web/packages/locfit/index.html      

## A. WGCNA Network Analysis ##  
### 1. Filtering and normalization of linear, effective transcript counts ###
Use the script `01_edgeR_normalization.R`. Inputs for this step are the two input files listed above (see *Required inputs to the repo include:* above).          

In brief, this script does the following:    
1) Imports data, sets up DGElist 
2) Filters data for low expression (i.e. requires X reads mapping in at least Y samples)
3) Normalizes data using TMM normalization
4) Estimates dispersions of counts
5) Produces normalized, log2 data as object `02_input_data/sfon_wgcna_01_output.RData` to be used in subsequent steps.
6) Produces mds plot of all samples


### 2. Set up data for WGCNA analysis, including subsets of samples and prepared traits ###
Use the script `02_WGCNA_setup.R`. Inputs for this step include the object `02_input_data/sfon_wgcna_01_output.RData` from the previous step. Much of this analysis comes from WGCNA package tutorials, although a lot is also custom or new code. Detailed instructions are found within the script.      

In brief, this script does the following:
1) Imports the data and the interpretation file, makes modifications to the traits and adds information for the Arctic Charr samples.   
2) Creates subsets of Brook Charr female offspring data, Brook Charr male offspring data, Arctic Charr male data, all sample backup
3) Produces all of these subsets into `02_input_data/sfon_wgcna_setup.RData`   

From here, either proceed to step 3 below (for full WGCNA analysis) or step 3A below, which is specific to plotting all Brook Charr transcriptome data with all selected trait data without using subsets.    

### 3. Build networks for female or male samples, and compare the opposite sex and Arctic Charr data to the network ###
Use the script `03_WGCNA_analysis.R`. Inputs for this step include the object `02_input_data/sfon_wgcna_setup.Rdata` from the previous step.  

To begin, one must choose what is the reference set to use to build the network, and the second set to compare to the reference network. Choose among "female", "male" and "AC" for female Brook Charr, male Brook Charr and Arctic Charr, respectively. Any combination of "female" or "male" as the network generation with any other as the module preservation test can be conducted. Here networks are not generated using Arctic Charr data.        

In brief, this script does the following: 
1) Filters based on low expression using only the samples within the reference set. (currently uses cpm >= 0.5, which is 14 transcripts mapping in the average library size (i.e. 27 M reads))    
2) Incorporates trait data, selects only the trait data required
3) Renames trait data and samples
4) Clusters samples based on expression and removes outlier samples
5) Calculates adjacency on all expressed transcripts, selects most connected 25 k transcripts, builds topological overlap matrix, generates module eigengenes, clusters module eigengenes and merges similar modules
6) Correlates module eigengenes with traits
7) Calculates module membership and gene significance for genes and traits of interest
8) Adds additional annotation information for all transcripts, save out (i.e. geneInfo)
9) Bring in second dataset, filter second dataset to be compared to network for low expression in the second dataset samples, keep only the genes that are present in the network
10) Test for module preservation of the network in the second dataset, and visualize results


#### Plot heatmaps of module genes
Use the script `01_scripts/plot_module_genes_in_heatmap.R` to plot heatmaps including genes from only one (or two) modules ordered by similarity in expression.  
Required inputs include:    
`04_results/<sex>_geneInfo_25k.csv` (which genes belong to which cluster in each sex)     
`03_normalized_data/normalized_output_matrix_log2.csv` (norm expr vals for all samples and genes)    
`00_archive/sfeq_interpretation_v1.csv` (sample metadata info)    


## B. Module and chromosome analysis ##  
### 1. Identify a single transcript per gene in the *de novo* reference transcriptome using a related species reference genome ###
In general terms, this maps a reference transcriptome against a reference genome and selects a single transcript for each contiguously overlapping alignment in the reference genome. Note: this has only been tested on the current dataset and has no guarantees for broader usage.   

*Required inputs to this section include:*
* reference genome (fasta)
* reference transcriptome (fasta) 

*Required software*
`gmap`    
`samtools`    
`bedtools`    

Use the following steps:   
##### i. Index reference genome with gmap #####      
For example, here with the Atlantic salmon genome (#cite)    
```
gmap_build -d ICSASG_v2 -D /home/ben/Documents/z-ssal_genome/ GCF_000233375.1_ICSASG_v2_genomic.fna
```

##### ii. Align transcripts against genome with gmap #####  
For example, here with the Phylofish Brook Charr reference transcriptome (#cite)   
```
gmap -D /home/ben/Documents/z-ssal_genome/ICSASG_v2 -d ICSASG_v2 -f samse -n 0 -t 4 ./sfontinalis_contigs_unwrap.fasta > ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam 2> ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam.log
```

##### iii. Filter mappings by quality and convert to bam #####   
```
samtools view -bS -q 30 ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam > ./sfontinalis_contigs_unwrap_v_ICSASG_v2_q30.bam     
samtools sort sfontinalis_contigs_unwrap_v_ICSASG_v2_q30.bam -o sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bam
samtools index sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bam`
```

##### iv. Convert bam alignments to bed file #####       
`bedtools bamtobed [OPTIONS] -i your_file.bam`

##### v. Determine sequence lengths for all transcripts #####      
```
cat sfontinalis_contigs_unwrap.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >  sfontinalis_contigs_unwrap_seq_lengths.txt
```
This approach was taken from http://www.danielecook.com/generate-fasta-sequence-lengths/      

##### vi. Identify a single transcript per contiguously mapping segment on reference genome #####  
Use the script `01_scripts/id_genes_from_ref_txome_w_ref_genome.R` to identify a single transcript per continuous mapping transcript segment.    
This script will require the following:    
* bed file of alignments (from step iv. above)
* transcript lengths file (from step v. above)
* sex-specific geneInfo file with added annotations, entitled <sex>_geneInfo_added_annot.txt (from A.3. above)  

In brief, this will:     
1. Link transcripts into contiguous overlapping segments, each of which will be considered one 'gene'. This produces <sex>_non-overlapping_transcripts.txt.       
2. Merge transcript length file with the non-overlapping transcripts file   
3. Merge above with the sex-specific geneInfo file with information about clusters to which each transcript belongs     
4. Select a single transcript to retain per 'unique gene' identifier from the non-overlapping transcripts. This ranks selection based on a) whether the transcript is expressed in the sex; b) the longest transcript. 

This outputs `<sex>_single_transcript_per_gene.txt`, which is the input for the next step.             

### 2. For each module, determine proportion of transcripts from each chromosome ###
In general terms, this characterizes the proportions of transcripts within each module that comes from each chromosome in the reference genome to see if there is any overrepresentation of specific chromosomes within specific modules.    
Use the script `plot_mod_chr_comp.R`, followed by `fishers_exact_test.R`      

*Required inputs to this section include:*
* sex-specific single transcript per contiguous alignment, from above '<sex>_single_transcript_per_gene.txt'. This file contains all details needed for this section. 
         
In brief, `plot_mod_chr_comp.R` will:
1. Identify which scaffolds are the main chromosomes for the genome (currently specific to Atlantic Salmon), generates the file 'chr_of_interest_list.txt'    
2. Determine how many genes are in each chromosome in all expressed genes (i.e. baseline = all modules, including grey and those not in the 25 k most connected genes (named as 'low.corr' here)    
3. Determine how many genes are in each chromosome in each module. Plot these in pie charts. 
This will result in the following file:  `<sex>_count_of_genes_per_chr_per_module.txt`. Use this file as an input to the next script below. 

In brief, `fishers_exact_test.R` will:
1. For each module, compare the proportion of genes in a chromosome to the proportion of baseline genes in this chromosome
#todo add more detail to this, and comment the code. 
This will result in the following file: `<sex>_fisher_test_chr_enrich.txt`   

More reading:    
http://www.biostathandbook.com/fishers.html    
http://www.dummies.com/education/math/statistics/how-to-compare-two-population-proportions/       
https://onlinecourses.science.psu.edu/stat414/node/268     


### 3. Plot positions of genes in enriched chromosomes ###
In general terms, this plots the positions of genes in the enriched module-chromosome combinations.   
*Required inputs to this section include:*
* output from B.2. above, <sex>_count_of_genes_per_chr_per_module.txt
* also from B.2. above, the list of the 29 chromosomes of Atlantic salmon: chr_of_interest_list.txt   
* the unwrapped reference genome used to align transcripts against above (see ## todo: link ## for unwrapping fasta)


Use the following steps:   
##### i. Determine full chromosome lengths for all chromosomes #####
Using the chr_of_interest_list.txt, determine chromosome lengths 
```
# select only chr of interest from unwrapped fasta:      
for i in $(cat chr_of_interest_list.txt) ; do grep -A 1 $i ~/Documents/z-ssal_genome/GCF_000233375.1_ICSASG_v2_genomic_unwrapped.fna ; done > chr_of_interest_genome.fasta     
# calculate chromosome lengths as done above for transcripts:   
cat chr_of_interest_genome.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > chr_of_interest_genome_lengths.txt        
```
This produces the file `chr_of_interest_genome_lengths.txt`

##### ii. Plot genes on chromosomes #####
Use the R script `plot_genes_on_chr.R` to generate plots.

This step requires:   
* Chromosome length file from above 'chr_of_interest_genome_lengths.txt' 
* Sex-specific transcriptome to genome alignment file '<sex>_single_transcript_per_gene.txt'
* Output of Fisher exact tests above '<sex>_fisher_test_chr_enrich.txt'

In brief, this plots the position of genes from the module of interest along the length of the chromosome of interest. This step chooses module-chromosome combinations based on the significant results from above.   





