library(ggplot2)
library(gridExtra)
library(ggpubr)

cbPalette2 <- c("#0072B2", "#D55E00", "#CC79A7","#56B4E9", "#009E73", "#F0E442")

# input 
# cluster members cds
# cluster member translated 

cluster_members_nucl

cluster_members_prot

groups



# Function to merge cluster finder output for protein and dna sequences.
merge_cluster_finder<-function(allele_clust, protein_clust){
  
  allele_info<-read.csv(allele_clust, header = FALSE) 
  protein_info<-read.csv(protein_clust, header = FALSE)
  
  # order both info to merge #
  df<-allele_info[-1,]
  df<-df[order(df$V2),]
  
  df2<-protein_info[-1,]
  df2<-df2[order(df2$V2),]
  
  # merge dfs #
  new<-cbind(df[1], df2)
  colnames(new)<-c("allele", "protein", "individual")
  
  return(new)
}


# Function to get Relative frequency (percentage of population with given genotype or phenotype)
out_table<-NULL
Count_perc<-function(x, prot_table, freq){
  attach(prot_table)
  total_n<-nrow(prot_table)
  for (i in seq(1:x)){
    n<-length(which(eval(freq) == i))
    count<-n[[1]]
    perc<-(count/total_n)*100
    ID<-i[[1]]
    Q<-cbind(ID, count, perc)
    out_table<-rbind(out_table, Q)
  }
  
  out_table<-data.frame(out_table)
  return(out_table)
}


Merged<-merge_cluster_finder(cluster_members_nucl, cluster_members_prot)





Host<-c(rep("D1", 10), rep("D2", 13), rep("D3", 15), rep("D4", 11), rep("D5", 8))

Host_allele_prot_table<-cbind(Allele_protein_members, Host)


Host_counts<-as.data.frame(table(Host_allele_prot_table$Host))
colnames(Host_counts)<- c("Host", "Freq")

# 
# # How many individuals with alleles (1:5) per host ... 
# Allele_frequency_all<-as.data.frame(table(Allele_protein_members$allele))
# 
# # For each host find percentage of individuals with each allele and make a variable which keeps
# # host info. 
# 
# XX<-NULL
# for(i in Host_counts$Host){
#   
#   # split into host groups
#   X<-Host_allele_prot_table[Host_allele_prot_table$Host == i,]
#   
#   # Get host ids for df - n ids should = n alleles
#   J<-data.frame(rep(i, 1))
#   colnames(J)<-"host"
#   
#   # Get percentage alleles
#   Q<-Count_perc(1, X, allele)   #allele is option (df colname)
#   S<-cbind(Q,J)
#   XX<-rbind(XX,S)
# }
# 
# All_host_percentage<-XX
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Notes 
# 
# # Full Method 
# #Aim: To plot Relative genotypic frequency and phenotypic frequency of host populations
# # Chromatograms and (Basecap)
# # Final MF file..all individuals.
# # CD63_all_individuals.fasta.
# # Translate CD63_all_individuals.fasta in Aliview, save as translated.
# # Run allele_finder on CD63_all_individuals.translated.fasta to get protein clusters.
# # sort allele and protein clusters by Individual ID.
# # Merge to get dataframe of both.
# # Plot Relative genotypic frequency for each host population and all populations.
# 
# 
# 
# 
# 
# 
# 
# 











