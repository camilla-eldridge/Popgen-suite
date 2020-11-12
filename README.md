Popgen suite

<br /> <br /> <br />

A suite of scripts for basic population genetic analysis from fasta sequences.

<br /> <br /> <br />

1. From a population sample of fasta sequences predict how many alleles are present and allele counts.

allele_finder

<br /> <br /> <br />

allele_finder.py


    Given multi-fasta file, sorts sequences into identical groups (can be used to cluster any strings in fasta format by 100% id)

    Returns tables of counts and cluster ids for plotting and further analysis.


       USAGE: allele_finder.py mf_file id


       Example: allele_finder.py test.fasta test
       
       
<br /> <br /> <br />

2. As allele_finder works on any string this tool can be used to predict protein variants and get protein counts in a population.











