library(WGCNA)
allowWGCNAThreads(nThreads=8)
options(stringsAsFactors=FALSE)


#set up data frame from expression data
load("/home/FCAM/tyankee/ANALYSIS/RNA-seq/Cranio_facial/CS/DESeq2/sva_norm_lg_cts_CF_FB_NCC.Rdata")



sva_norm_lg_cts_t<-t(sva_norm_lg_cts)

rownames(sva_norm_lg_cts_t)<-colnames(sva_norm_lg_cts)

#declare variables
npower<-18

#Make TOM
adjacency=adjacency(sva_norm_lg_cts_t,power=npower,type="signed")
TOM=TOMsimilarity(adjacency,TOMType="signed")
dissTOM=1-TOM

save(dissTOM,file="dissTOM.Rdata")

#cluster using TOM 
geneTree=hclust(as.dist(dissTOM),method="average")

save(geneTree,file="gene_Tree.Rdata")
