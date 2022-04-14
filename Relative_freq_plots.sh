#!/bin/bash

#author: camilla eldridge

nuc_seq="$1"
prot_seq="$2"
id="$3"
groups="$4"
group_id="$5"

# Notes on usage #

if [ "$1" == "-h" ]; then
  echo "Usage:  ./Relative_freq_plot.sh  nucl.fasta prot.fasta gene_id  groups.csv group_id"
  exit 0
fi

# Get allele count and cluster info from population sequences #
# (each sequence represents and individual #
cluster_finder.py "$nuc_seq" "$id"_nucl

# Get protein countsand cluster info (prot_seq = translated nuc_seq sequences) #
cluster_finder.py "$prot_seq" "$id"_prot

# Run script to plot Relative frequencies
./Relative_geno_pheno_plots.R "$id"_nucl_cluster_members.csv "$id"_prot_cluster_members.csv "$id"_nucl_cluster_info.txt "$id"_prot_cluster_info.txt "$groups" "$group_id" "$nuc_seq"


