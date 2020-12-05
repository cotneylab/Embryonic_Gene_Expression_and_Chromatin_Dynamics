##Tutorial for clusterProfiler package in R 

library(clusterProfiler)
library(org.Hs.eg.db)


##Convert gene ids to entrez 
##If your ensembl IDs have ".XX" suffix, you need to remove using this code: 
#gsub("\\..*","",df$geneID)

##genes to GO
cl1_up<-mapIds(org.Hs.eg.db,keys=genes2$cl2,column="ENTREZID",keytype="ENSEMBL",multiVals="first")

##background set of genes 
uni<-mapIds(org.Hs.eg.db,keys=genes2$cl2,column="ENTREZID",keytype="ENSEMBL",multiVals="first")


##GO over-representation test 

ego = enrichGO(gene=cl1_up$ENTREZD, universe=uni$ENTREZID,
	OrgDb=org.Hs.eg.db, ont ="BP",
	pAdjustMethod= "BH",
	pvalueCutoff =0.01,
	qvalueCutoff=0.05,
	readable=TRUE)

#dotplot 

