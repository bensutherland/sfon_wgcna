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

## 3. Identify one isoform per gene in reference transcriptome using a reference genome ##
First, index your reference genome for use with gmap, for example, here with S. salar.       
`gmap_build -d ICSASG_v2 -D /home/ben/Documents/z-ssal_genome/ GCF_000233375.1_ICSASG_v2_genomic.fna`

Second, align against the reference genome
`gmap -D /home/ben/Documents/z-ssal_genome/ICSASG_v2 -d ICSASG_v2 -f samse -n 0 -t 4 ./sfontinalis_contigs_unwrap.fasta > ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam 2> ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam.log`

Filter based on quality and convert to bam    
```
samtools view -bS -q 30 ./sfontinalis_contigs_unwrap_v_ICSASG_v2.sam > ./sfontinalis_contigs_unwrap_v_ICSASG_v2_q30.bam     
samtools sort sfontinalis_contigs_unwrap_v_ICSASG_v2_q30.bam -o sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bam
samtools index sfontinalis_contigs_unwrap_v_ICSASG_v2_q30_sorted.bam`
```

Then finally create a bed file     
`bedtools bamtobed [OPTIONS] -i your_file.bam`

Note: in case you need to find sequence lengths:    
`cat sfontinalis_contigs_unwrap.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' >  sfontinalis_contigs_unwrap_seq_lengths.txt`
This was from: http://www.danielecook.com/generate-fasta-sequence-lengths/

This bed file will be used as an input to the following R script to identify a single transcript per putative gene:     
`id_genes_from_ref_txome_w_ref_genome.R`     
In this R script, you will:   
1. Link transcripts into contiguous overlapping segments, each of which will be considered one 'gene'    
2. Associate fasta accession lengths to the 'unique gene' identifier.  
3. Import a sex-specific geneInfo object that contains information about the clusters to which each gene belongs.     
4. Select a single transcript to retain per 'unique gene' identifier in order of a) is expressed and; b) take longest. Possible to export list of single transcripts with their clusters here.   
5. Plot each clusters' proportions of each Atlantic Salmon chromosome.    

## 4. Identify paralogs 
A reciprocal best hit blast has been constructed to blast same-on-same but to remove the focal transcript from the blast so that the first transcript can be chosen as the 'best hit' without hitting itself. 

This functions by using BLAST of focal transcript against all reference transcriptome except itself, taking the single top hit and writing onto 'blast_output.txt'     

Use the following script to perform within transcriptome RBH blast (removing focal transcript and use BLAST against itself:    
`repeated_blast.sh`    

First make a blast database with the unwrapped reference transcriptome:   
`makeblastdb -in ./GCF_000233375.1_ICSASG_v2_genomic.fna -dbtype nucl -parse_seqid`   

Input: unwrapped ref transcriptome fasta, file with 'all_accn_names.txt'       
Produces 'blast_output.txt'   

Next need to use 'runall.sh' but re-code that does not assume the first best hit would be to itself, as these have been removed.    


Once this is constructed, the paralog pairs identified can be associated with the object from 3.4 export above.    






Next, a script must be written to identify whether paralogs are in the same or different cluster, then plotted number same vs number different and tested statistically if necessary.     

