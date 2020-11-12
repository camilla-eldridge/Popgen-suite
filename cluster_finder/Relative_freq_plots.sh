#!/bin/bash

nuc_seq="$1"
prot_seq="$2"
id="$3"

# Notes on usage #

if [ "$1" == "-h" ]; then
  echo "Usage:  ./Relative_freq_plot.sh  nucl.fasta prot.fasta gene_id  translation_table"
  exit 0
fi

# Get allele count and cluster info from population sequences #
# (each sequence represents and individual #
cluster_finder.py "$nuc_seq" "$id"_nucl

# Get protein countsand cluster info (prot_seq = translated nuc_seq sequences) #
cluster_finder.py "$prot_seq" "$id"_prot



