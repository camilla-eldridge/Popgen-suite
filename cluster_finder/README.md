# Cluster finder 

* Given multi-fasta file, sorts sequences into identical groups (can be used to cluster any strings in fasta format by 100% id)

* Returns tables of counts and cluster ids for plotting and further analysis. <br /> <br /> <br /> 



         USAGE: cluster_finder.py mf_file id


         Example: cluster_finder.py test.fasta test


 <br /> <br /> <br /> 
 
Output: <br /> <br /> <br /> 

         test_clusters.fasta  

Sequence file sorted into indentical cluster groups, n in cluster shown in header.  <br /> <br /> <br /> 

         test_cluster_info.txt  
Table of cluster members, n in cluster and cluster representative id    <br /> <br /> <br /> 

         test_cluster_members.csv  

Table of cluster ids per individual   <br /> <br /> <br /> 

         test_cluster_reps.fasta
         
Cluster representative sequences, in order (cluster 1,2,3...)
