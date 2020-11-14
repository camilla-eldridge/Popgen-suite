#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(gridExtra)
library(ggpubr)

# Input as csv output from cluster_finder.py #
cluster_members_nucl<-args[1]
cluster_members_prot<-args[2]
cluster_info_nucl<-args[3]
cluster_info_prot<-args[4]
groups<-args[5]
group_name<-args[6]
# 
# # Allele ID for each individual in population.
# cluster_members_prot<-"test_nucl_cluster_members.csv"
# 
# # Protein ID for each individual in population.
# cluster_members_nucl<-"test_prot_cluster_members.csv"
# 
# # Whole pop allele counts.
# cluster_info_nucl<-"test_nucl_cluster_info.txt" 
# 
# # Whole pop protein counts. 
# cluster_info_prot<-"test_prot_cluster_info.txt"
# 
# # Info on group assignment for each individual 
# groups<-"group_info.csv"
# 
# # name of group for plotting and table generation
# group_name="Host"


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

attach(reord_A
       )
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


# arrange in grid 
ggarrange(Relative_allele_plot,Relative_prot_plot, labels = c('I','II'), hjust = 0.1, vjust =1.2, ncol = 2, nrow = 1)


# save plots
ggsave("Relative_phen_geno_plot.png", device = "png", width = 25, height = 10, units = "cm")


