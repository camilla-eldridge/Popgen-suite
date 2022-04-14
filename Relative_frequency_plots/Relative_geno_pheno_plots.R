#!/usr/bin/env Rscript

#author: camilla eldridge
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(pegas)
require(ape)


# Input as csv output from cluster_finder.py #
cluster_members_nucl<-args[1]
cluster_members_prot<-args[2]
cluster_info_nucl<-args[3]
cluster_info_prot<-args[4]
groups<-args[5]
group_name<-args[6]
cds_seq<-args[7]

# # Allele ID for each individual in population.
#cluster_members_prot<-"test_nucl_cluster_members.csv"
# 
# # Protein ID for each individual in population.
#cluster_members_nucl<-"test_prot_cluster_members.csv"
# 
# # Whole pop allele counts.
#cluster_info_nucl<-"test_nucl_cluster_info.txt" 
# 
# # Whole pop protein counts. 
#cluster_info_prot<-"test_prot_cluster_info.txt"
# 
# # Info on group assignment for each individual 
#groups<-"group_info.csv"
# 
# # name of group for plotting and table generation
#group_name="Host"

# fasta sequences  
#cds_seq<-"test.fasta"


####################        Part I          #################################

# To output a table (csv) displaying allele, protein 
# and group infor for each sequence in the population sample.

# Function to merge cluster finder output for protein and dna sequences.
merge_cluster_finder<-function(allele_clust, protein_clust){
  
  allele_info<-read.csv(allele_clust, header = FALSE) 
  protein_info<-read.csv(protein_clust, header = FALSE)
  
  # order by sequence id to merge #
  df<-allele_info[-1,]
  df<-df[order(df$V2),]
  
  df2<-protein_info[-1,]
  df2<-df2[order(df2$V2),]
  
  # merge dfs #
  new<-cbind(df[1], df2)
  colnames(new)<-c("allele", "protein", "individual")
  
  return(new)
}

# Merge protein and allele info for each individual #
Merged<-data.frame(merge_cluster_finder(cluster_members_nucl, cluster_members_prot))

# read in group allocations #
Seq_groups<-data.frame(read.csv(groups, header = TRUE))

# order table by Individual to merge group if with Merged table of protein and allele ID's #
Ordered_groups<-Seq_groups[order(Seq_groups$individual),]

# Merge tables #
Merged_ordered<-cbind(Merged, Ordered_groups$group)

# rename colomns, adding user info for group name # 
colnames(Merged_ordered)<-c("Allele", "Protein", "Individual", paste(group_name))

# write the final merged table to file # 
result <- write.csv(x=Merged_ordered, file="Allele_and_protein_variants.csv", row.names = FALSE)

# get total n alleles and protein variants in population sample #
n_alleles<-as.numeric(length(unique(Merged_ordered$Allele)))
n_proteins<-as.numeric(length(unique(Merged_ordered$Protein)))

# get unique groups  #
groups_uniq<-unique(Seq_groups$group)

# get allele and protein variant counts for entire population sample # 
all_allele_counts<-read.csv(cluster_info_nucl, header = TRUE)
all_protein_counts<-read.csv(cluster_info_prot, header = TRUE)

# Get number of individuals in sample #
N_individuals<-sum(all_allele_counts$n_in_cluster)

# Percentage function #
Percentage<-function(N, counts){
  all_perc=NULL
  X<-split(counts, sort(as.numeric(rownames(counts)))) # splits colomns by row
  for (i in X){             # for each colomn in each row... #
    perc<-(i[[2]]/N)*100
    count_perc<-cbind(i[[1]], i[[2]], perc)
    all_perc<-rbind(all_perc, count_perc)
  }
  return(all_perc)
}


# Make group labels colomn for all alleles and prots #
group_all_alleles<-data.frame(rep("All", n_alleles))
group_all_prots<-data.frame(rep("All", n_proteins))

# Get percentage of alleles in entire population #
Perc_alleles<-cbind((Percentage(N_individuals, all_allele_counts)), group_all_alleles)
colnames(Perc_alleles)<-c("Allele", "Count", "Perc_alleles", paste(group_name))

# Get percentage of protein variants in entire population #
Perc_proteins<-cbind((Percentage(N_individuals, all_protein_counts)), group_all_prots)
colnames(Perc_proteins)<-c("Protein", "Count", "Perc_proteins", paste(group_name))


####################        Part II          #################################

# Get percentage of alleles and protein variants in each group 
# Get number of individuals in sample #
# N_individuals<-sum(all_allele_counts$n_in_cluster)


# for prot and allele percs
Allele_group_percs<-NULL
Protein_group_percs<-NULL


for(j in groups_uniq){
  f<-paste(group_name)
  X2<-Merged_ordered[Merged_ordered[[f]] == j,] 
  
  # n individuals in group 
  N_ind_group<-length((X2$Individual))
  
  # Get group name 
  gr_name<-X2[4][1,]
  
  # Get counts for allele and prot variants
  count_allele<-as.data.frame(table(X2$Allele, exclude = "Cluster ID")) # stop cluster id from being a level.
  colnames(count_allele)<-c("Allele", "count")
  
  count_prot<-as.data.frame(table(X2$Protein, exclude = "Cluster ID"))
  colnames(count_prot)<-c("Protein", "count")
  
  # Get percs for allele and prot variants
  group_perc_allele<-Percentage(N_ind_group, count_allele)
  group_perc_prot<-Percentage(N_ind_group, count_prot)
  
  # add colomn of host id
  host_gr_a<-data.frame(rep(as.factor(gr_name), n_alleles))
  group_allele_table<-cbind(host_gr_a, group_perc_allele)
  
  host_gr_p<-data.frame(rep(as.factor(gr_name), n_alleles))
  group_prot_table<-cbind(host_gr_p, group_perc_prot)

  # edit col names #
  colnames(group_allele_table)<-c(paste(group_name), "Allele", "Count", "Perc_alleles")
  colnames(group_prot_table)<-c(paste(group_name), "Protein", "Count", "Perc_proteins" )
  
  Allele_group_percs<-rbind(Allele_group_percs, group_allele_table)
  Protein_group_percs<-rbind(Protein_group_percs, group_prot_table)
  
  }


#  write out per group counts and percs 
result2 <- write.csv(x=Allele_group_percs, file="group_allele_counts.csv", row.names = FALSE)
result3 <- write.csv(x=Protein_group_percs, file="group_protein_counts.csv", row.names = FALSE)



# ####################        Part III          #################################
# colour blind palette for gglot2
cbPalette <- c("#0072B2", "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442")
 
# reorder by perc for plotting
reord_A<-Allele_group_percs[order(Allele_group_percs$Perc_alleles),]

attach(reord_A)

# string for x axis variable 
x_axis<- paste(group_name)

# make plot for alleles
Relative_allele_plot<-ggplot() +
    geom_bar(data = reord_A, aes_string(x = x_axis, y = reord_A$Perc_alleles, fill = as.factor(Allele)), stat = "identity") +
  theme_classic() + ylab("Relative genotypic frequency (%)") + xlab(paste(group_name)) +
  theme(axis.text.x = element_text(size = 12)) + theme(axis.title = element_text(size = 15)) + theme(axis.text.y = element_text(size =12), axis.title.y=element_text(size=15)) +
  scale_fill_manual("Allele", values = cbPalette)

detach(reord_A)


# make plot for proteins
# reorder by perc for plotting
reord_P<-Protein_group_percs[order(Protein_group_percs$Perc_proteins),]

attach(reord_P)

Relative_prot_plot<-ggplot() +
  geom_bar(data = reord_P, aes_string(x = x_axis, y = reord_P$Perc_proteins, fill = as.factor(Protein)), stat = "identity") +
  theme_classic() + ylab("Relative phenotypic frequency (%)") + xlab(paste(group_name)) +
  theme(axis.text.x = element_text(size = 12)) + theme(axis.title = element_text(size = 15)) + theme(axis.text.y = element_text(size =12), axis.title.y=element_text(size=15)) +
  scale_fill_manual("Protein", values = cbPalette)




########### AMOVA, TAJD ############
library(pegas)

# All individuals sequenced #
indiv_seqs<-read.dna(cds_seq, format = "fasta", skip = 0, as.character = FALSE, as.matrix = TRUE)
class(indiv_seqs)

#Tajima with p value (pegas) 
t<-tajima.test(indiv_seqs)

# gc content
gc<-GC.content(indiv_seqs)

#nucleotide diversity
d<-nuc.div(indiv_seqs, variance = FALSE, pairwise.deletion = FALSE)

#haplotype diversity
h<-hap.div(indiv_seqs)

#haplotype number
hn<-haplotype(indiv_seqs)

# ape to find base frequencies #
bf<-base.freq(indiv_seqs)

# ramos onsins rozas test
Ra<-R2.test(indiv_seqs)



# plot arrangement
jpeg("summary.jpg")

par(mfrow = c(2,2))

# call r2 plot 
#Ra
R2.test(indiv_seqs)

# plot of haplotype freqs
plot(haplotype(indiv_seqs), main = "Haplotype distribution")

# Site frequency spectrum
sp<-site.spectrum(indiv_seqs)
plot(sp, col = "lightgreen", main = "SFS")
#Wakeley, J. (2009) Coalescent Theory: An Introduction. Greenwood Village, CO: Roberts and Company Publishers.

# R2 test -  Ramos-Onsins-Rozas Test of Neutrality
#R2.test(indiv_seqs)
#Ramos-Onsins, R. and Rozas, R. (2002) Statistical properties of new neutrality tests against population growth. Molecular Biology and Evolution, 19, 2092–2100.
#Sano, J. and Tachida, G. (2005) Gene genealogy and properties of test statistics of neutrality under population growth. Genetics, 169, 1687–1697.

# pairwise mismatch 
MMD(indiv_seqs, col = "grey", breaks = 20, legend = FALSE, lty = 1:2, main = "Pairwise mismatch")
# blue = empirical
# red = stable expectation
#Rogers, A. R. and Harpending, H. (1992) Population growth makes waves in the distribution of pairwise genetic-differences. Molecular Biology and Evolution, 9, 552–569.

dev.off()

######### summary tables ######

df.summary<-cbind(t[1], t[2],d,h,Ra[1], Ra[2],bf[1],bf[2], bf[3], bf[4])
colnames(df.summary)<-c("D", "D.pval", "nuc.div", "hap.div", "r2", "r2.pval", "a", "c", "g", "t")
row.names(df.summary)<-NULL

############  AMOVA #############################
## genetic distance matrix for all individuals ##
indiv_seqs.dist <- dist(x = indiv_seqs, method = "euclidean", diag=TRUE, upper=TRUE)

# order of input sequences stays the same so we can use Host variable for AMOVA.
Host<-Merged_ordered$Host  # from groups<-"group_info.csv"

AMOVA<-amova(indiv_seqs.dist ~ Host, nperm = 7000)

#sum of squares deviation, mean square deviation 
AMOVA$tab

# variancecoeff (a)
AMOVA$varcoef

# variance componants (sigma2, pvalue)
AMOVA$varcomp

#phiST
PhiST <- 1 - AMOVA$varcoef/(AMOVA$varcoef+AMOVA$varcomp[2,1])


# AMOVA tables for ggplot #
df.amova<-data.frame(AMOVA$varcomp)
ss <- tableGrob(df.amova)

f<-tableGrob(AMOVA$tab)

phi.varco<-rbind(as.numeric(PhiST), as.numeric(AMOVA$varcoef))
rownames(phi.varco)<-c("PhiST", "var.coef")
colnames(phi.varco)<-NULL

k <- tableGrob(phi.varco)

sums<-tableGrob(df.summary)



# arrange in grid 
p1<-ggarrange(Relative_allele_plot,Relative_prot_plot, labels = c('I','II'), hjust = 0.1, vjust =1.2, ncol = 2, nrow = 1)
p2<-ggarrange(ss,f,k, ncol = 1, nrow = 3, labels = "III")
p3<-

ggarrange(p1, p2, ncol = 1, nrow = 2)


# save plots
ggsave("Relative_phen_geno_plot_amova.png", device = "png", width = 25, height = 25, units = "cm")


