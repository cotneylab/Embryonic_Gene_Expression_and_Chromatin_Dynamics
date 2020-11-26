
diff_gene_list<-vector(mode="list", length=4)
names(diff_gene_list)<-c("CS13","CS14","CS15","CS17")

cs13<-c("CS13CS14","CS13CS15","CS13CS17")
cs14<-c("CS14CS13","CS14CS15","CS14CS17")
cs15<-c("CS15CS13","CS15CS14","CS15CS17")
cs17<-c("CS17CS13","CS17CS14","CS17CS15")

li_DE<-list("cs13"=cs13,"cs14"=cs14,"cs15"=cs15,"cs17"=cs17)

for(j in c(1:4)){
df_tmplist<-NULL
  
for(i in li_DE[[j]]){

df<-read.csv(paste0(i,"_DE.tsv"),sep="\t")
df_sub<-subset(df, log2FoldChange > 1 | log2FoldChange < -1 & padj < 0.05)
df_tmplist<-c(df_tmplist,row.names(df_sub))

}
diff_gene_list[[j]]<-df_tmplist
}