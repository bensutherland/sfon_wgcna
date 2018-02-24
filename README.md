# **WGCNA Analysis for Brook Charr** #
WGCNA analysis of transcriptome data for Sutherland et al. 2018 (in prep).   
This pipeline has only been tested and used for the purposes of the authors of the above article, and has no guarantees of usefulness whatsoever for other purposes.   

*Required inputs to the repo include:*   
1) a counts file with linear values that are output from the following repo:  
https://github.com/bensutherland/SE-reads_assemble-to-counts.git    
specifically using `01_script/03_express.sh` and then `01_scripts/utility_scripts/prepare_gxlevels_matrix.R`.   
Essentially, this uses the 'effective counts' column as suggested by eXpress (#cite).
This should be entitled/located as follows: `02_input_data/out.matrix.csv`

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

These steps are performed in R using the WGCNA package. See detailed instructions within `01_scripts/03_WGCNA_BC_fem_mods_male_comp.R`     


## 4. Identify one isoform per gene in reference transcriptome using a reference genome ##
This section will allow you to, given a reference transcriptome and genome input, identify which transcripts map and overlap in position on the reference genome. Then it allows you to retain only one transcript per contiguous map mapping transcript unit. (no guarantees, experimental)          

Required: reference genome, reference transcriptome   

First, index reference genome gmap, e.g. with S. salar:       
`gmap_build -d ICSASG_v2 -D /home/ben/Documents/z-ssal_genome/ GCF_000233375.1_ICSASG_v2_genomic.fna`

Second, align transcripts to reference genome:
`gmap -D /home/ben/Documents/z-ssal_genome/ICSASG_v2 -d ICSASG_v2 -f samse -n 0 -t 4 ./sfontinalis_contigs_unwrap.fasta > ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam 2> ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam.log`

Third, filter mappings by mapq and convert to bam:    
```
samtools view -bS -q 30 ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam > ./sfontinalis_contigs_unwrap_v_ICSASG_v2_q30.bam     
samtools sort sfontinalis_contigs_unwrap_v_ICSASG_v2_q30.bam -o sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bam
samtools index sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bam`
```

Fourth, create a bed file     
`bedtools bamtobed [OPTIONS] -i your_file.bam`

Fifth, find sequence lengths of all transcripts:    
`cat sfontinalis_contigs_unwrap.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >  sfontinalis_contigs_unwrap_seq_lengths.txt`
This was from: http://www.danielecook.com/generate-fasta-sequence-lengths/

Sixth, use R script to identify a single transcript per continuous mapping transcript segment: `id_genes_from_ref_txome_w_ref_genome.R`.    
This Rscript will require: a bed file (from above), a transcript lengths file, a <sex>_geneInfo_added_annot.txt file.         
 The steps required are as follows:      
1. Link transcripts into contiguous overlapping segments, each of which will be considered one 'gene'    
2. Associate transcript lengths to the 'unique gene' identifier.  
3. Import a sex-specific geneInfo object (see above) that contains information about the clusters to which each gene belongs.     
4. Select a single transcript to retain per 'unique gene' identifier in order of a) is expressed and; b) take longest. Possible to export list of single transcripts with their clusters here.   

This data is output as a table, such as `<sex>_single_transcript_per_gene.txt`.     
Move to next step to plot modules' chromosomal composition.     


## 5. Compare representation by chromosome across modules ##
Use Rscript to characterize the proportions of genes from each chromosome in the modules, and compare to the baseline of all genes, using the following script: `plot_mod_chr_comp.R`         
Basically, this will determine which scaffolds are the chromosomes, what the baseline is in terms of number of genes from each chromosome, calculate this for each module, plot in pie charts and export to text the counts of genes from each chromosome in each module.    

#todo: Apply Fisher's exact test and bonferonni correction.   


## 6. Plot positions of genes in enriched chromosomes ##
This uses as input the file created from `plot_mod_chr_comp.R`, `chr_of_interest_list.txt`, and the reference genome (unwrapped).         

Determine chromosome lengths of 'chr of interest', first selecting only the chr of interest:   
`for i in $(cat chr_of_interest_list.txt) ; do grep -A 1 $i ~/Documents/z-ssal_genome/GCF_000233375.1_ICSASG_v2_genomic_unwrapped.fna ; done > chr_of_interest_genome.fasta`    

Then use the same command as above to calculate lengths of the accessions in this fasta file:    
`cat chr_of_interest_genome.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > chr_of_interest_genome_lengths.txt`    


Use Rscript to do plotting `plot_genes_on_chr.R`    







## 4. Identify paralogs 
A reciprocal best hit blast has been constructed to blast same-on-same but to remove the focal transcript from the blast so that the first transcript can be chosen as the 'best hit' without hitting itself. 

This functions by using BLAST of focal transcript against all reference transcriptome except itself, taking the single top hit and writing onto 'blast_output.txt'     

Use the following script to perform within transcriptome RBH blast (removing focal transcript and use BLAST against itself:    
`repeated_blast.sh`    

First make a blast database with the unwrapped reference transcriptome:   
`makeblastdb -in ./GCF_000233375.1_ICSASG_v2_genomic.fna -dbtype nucl -parse_seqid`   

Input: unwrapped ref transcriptome fasta, file with 'all_accn_names.txt'       
Produces 'blast_output.txt'   

Next need to use 'runall.sh' steps but re-code that does not assume the first best hit would be to itself, as these have been removed.    

This will probably require automation bash script to perform the necessary steps.    
Once this is constructed, the paralog pairs identified can be associated with the object from 3.4 export above.    

As a final step, a script must be written to identify whether paralogs are in the same or different cluster, then plotted number same vs number different and tested statistically if necessary.     

