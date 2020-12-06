##Tutorial for clusterProfiler package in R 

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(reshape2)
library(ggraph)
library(igraph)



#Also source all code from
#https://github.com/YuLab-SMU/enrichplot/blob/master/R/emapplot.R  (use the emapplot_reconfig.R that has edited code tara made)
#https://github.com/YuLab-SMU/enrichplot/blob/master/R/utilities.R

##If your ensembl IDs have ".XX" suffix, you need to remove using this code: 
#gsub("\\..*","",df$geneID)
#make list of genes per cluster
#HM_genes<-list("pink"=pink,"purple"=gsub("\\..*","",row.names(clus_genes$cl2_purple)))

##GO over-representation test 


enrich_GO_HM<-vector(mode="list",length=6)


#Input list for each cluster of genes

enrich_GO_HM<-function(genes){

  enrich_GO_HM<-vector(mode="list",length=3)

ontology<-c("MF","CC","BP")  

names(enrich_GO_HM)<-ontology

for( i in ontology){

  ego_tmp<-enrichGO(genes,keyType = "ENSEMBL",ont=i,OrgDb = org.Hs.eg.db)

ego_tmp@result$plog<--log(ego_tmp@result$p.adjust,2)

   enrich_GO_HM[[i]]<-ego_tmp

}

  return(enrich_GO_HM)

}

li_enrichGO<-lapply(HM_genes,enrich_GO_HM)
#To plot, change colors as needed in the emapplot_reconfig.R and source.
#source("emapplot_reconfig.R")
#emapplot(li_enrichGO$pink$BP)+theme_bw()



