##Tutorial for clusterProfiler package in R 

library(clusterProfiler)
library(org.Hs.eg.db)


##Convert gene ids to entrez 
##If your ensembl IDs have ".XX" suffix, you need to remove using this code: 
#gsub("\\..*","",df$geneID)

##genes to GO
cl1_up<-bitr(row.names(cl1), fromType='ENSEMBL', toType='ENTREZID', OrgDb = 'org.Hs.eg.db')

##background set of genes 
uni<-bitr(row.names(all), fromType='ENSEMBL', toType='ENTREZID', OrgDb = 'org.Hs.eg.db')


##GO over-representation test 

ego = enrichGO(gene=cl1_up$ENTREZD, universe=uni$ENTREZID,
	OrgDb=org.Hs.eg.db, ont ="BP",
	pAdjustMethod= "BH",
	pvalueCutoff =0.01,
	qvalueCutoff=0.05,
	readable=TRUE)
