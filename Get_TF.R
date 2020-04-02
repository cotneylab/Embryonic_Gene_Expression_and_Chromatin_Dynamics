#This function is used to extract TF from a list provided by david GO FuncAnnotChart outputs.
#It is currently set up to get TFs per cluster and output TFs by cluster.  (snatac data)

get_TF<-function(clusli){
df<-c() 
for(j in clusli){
  print(j)
    
cl_func<-read.table(paste0(j,"FuncAnnotChart.tsv"),header=TRUE,sep="\t")
cl_func_filt<-cl_func[grep("transcription factor",cl_func$Term),]
if(nrow(cl_func_filt) < 1){
  next
}

genes<-c()
for(i in c(1:nrow(cl_func_filt))){
  
  genes<-paste(genes, cl_func_filt$Genes[i])
  
}

genes<-strsplit(genes," ")
genes<-strsplit(genes[[1]],",")
genes<-list.rbind(genes)
genes<-as.data.frame(genes)
genes<-unique(genes$V1)

sym<-c()
for(k in genes){
  
  tmp<-geneGini_scatac_log_filt[which(geneGini_scatac_log_filt$ENS == k),]
  
  sym<-c(sym,as.character(tmp$Gene))
}

dftmp<-data.frame("TF"=sym,"ENST"=as.character(genes),"cl"=c(rep(j,length(sym))))
df<-rbind(df,dftmp)
}
  
  return(df)
}