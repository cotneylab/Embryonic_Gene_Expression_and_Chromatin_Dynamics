library("AnnotationDbi")
library("org.Hs.eg.db")

C<-c("NCC","CS13","CS14","CS15","CS17","CS22")


for(i in C){
print(paste("sample",i))
  CSsamp<-NULL
  UPDN<-NULL

  for(j in C){

    if (i==j){
	next
    }

    restmp<-results(ddssva,alpha=0.05,contrast=c("sample",i,j))

    summary(restmp)
    restmp$Symbol<-mapIds(org.Hs.eg.db,keys=gsub("\\..*","",row.names(restmp)),column="SYMBOL",keytype="ENSEMBL",multivals="first")
    restmp<-cbind(restmp,counts(ddssva))
    write.table(restmp,file=paste(i,j,"_DE.tsv",sep=""),sep="\t")

    CSsamp<-c(CSsamp,rep(j,2))

    UP<-subset(restmp,log2FoldChange >= 1 & padj <0.05)
    DN<-subset(restmp,log2FoldChange <= -1 & padj <0.05) 
    UPDN<-c(UPDN,nrow(UP),-nrow(DN))

  }

  Col<-c(rep(c("blue","red"),(length(C)-1)))

  pairwise<-data.frame("CSsamp"=CSsamp,"UPDN"=UPDN,"Col"=Col)
  colnames(pairwise)[2]<-paste(i,"_UPDN",sep="")
  write.csv(pairwise,file=paste(i,"vAll.csv",sep=""),row.names=FALSE)
}
