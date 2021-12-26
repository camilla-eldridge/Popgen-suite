# Popgen suite

A suite of scripts for basic population genetic analysis from fasta sequences that can be implemented in analysis pipelines.
<br /> <br /> <br />


Problem: Preparation of input files for popgen analyses can be laborious, particularly if you only want to run basic analyses, and most popgen software cannot be easily implimented into automated pipelines. 

Solution: Popgen-suite makes it easy to run basic population statistics, requiring a fasta file and a csv file to define population groups. Analyses can be carried out using separate scripts, enabling integration into automated pipelines.

<br /> <br /> <br />


## Allele and protein variant prediction.<br /> <br /> <br />

 
 1. Predicting Allele and Protein variants. <br /> <br /> <br />  

`cluster_finder.py` clusters sequences based on 100% identity and can be used on any kind of sequence or string. <br /> <br /> <br />

From a population sample of fasta sequences cluster_finder.py will predict how many alleles are present and allele counts. Running the same script 
on the translated sequences will predict the protein variants and their counts.<br /> <br /> <br />


<br /> <br /> <br />

`cluster_finder.py`


- Given multi-fasta file, sorts sequences into identical groups (can be used to cluster any strings in fasta format by 100% id).<br /> <br /> <br />
- Returns tables of counts and cluster ids for plotting and further analysis.<br /> <br /> <br />


       Example: cluster_finder.py   test.fasta   test
       
       
<br /> <br /> <br />
<br /> <br /> <br />


2. Plotting Relative Allele and Protein variant frequencies in multiple populations (e.g different host populations or geographic regions).

After running `cluster_finder.py` on both protein and nucleotide data sets, the output from both predictions can be merged.<br /> <br /> <br />
Running `Relative_geno_pheno_plots.R` will sort this output so that a comparative plot of Relative genotypic and phenotypic frequencies is produced.<br /> <br /> <br />

    Example of usage: Relative_geno_pheno_plots.R  test_nucl_cluster_members.csv  test_prot_cluster_members.csv   group_info.csv 
    
    
<br /> <br /> <br />
<br /> <br /> <br />
<br /> <br /> <br />


![alt text](Example/Relative_phen_geno_plot.png)


<br /> <br /> <br />

*To test for any evidence in population structure observed in the final plot an AMOVA can be calculated using pegas (see Population statistic summary tables).<br /> <br /> <br />


<br /> <br /> <br />
<br /> <br /> <br />

## Automated prediction and Relative frequency plotting. 

In `Relative_freq_plots.sh` steps 1 and 2 above are combined. 
<br /> <br /> <br />


     Example of usage: ./Relative_freq_plots.sh   test.fasta   test.translated.fasta   test   group_info.csv   Host
<br /> <br /> <br />

Notes on usage:
<br /> <br /> <br />
- Column names must be the same as shown in the test files (see group_info.csv). <br /> <br /> <br />
- Relative frequencies and sorted allele and protein variation information for each individual in the dataset is always ouput as:<br /> <br /> <br />
`group_allele_counts.csv`, `group_protein_counts.csv` and `Allele_and_protein_variants.csv` <br /> <br /> <br />

<br /> <br /> <br />
<br /> <br /> <br />


## Population statistic summary tables. 

coming soon...
<br /> <br /> <br />

NB. add in :

PhiST <- 1 - fit$varcoef/(fit$varcoef+fit$varcomp[2,1])
PhiST


## Per site analyses of Tajima's D and Pi.

coming soon...
<br /> <br /> <br />



