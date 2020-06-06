###

# The input is expected to be an adjacency matrix. The rows and columns are the same 
# The diagnol is equal to 1. 
#cutoff is a value between 0 and 1.
#col_mod is a string of the module color 
#The output is a dataframe for inputting into cytoscape network. 
library(AnnotationDbi)
library(org.Hs.eg.db)


makenetwork<-function(inputmat,cutoff,col_mod){
  
  finaldf<-data.frame("Source"="","Target"="","cor"="")
  finaldf<-finaldf[FALSE,]
  n<-1
  
  for (i in row.names(inputmat)){
    n<-n+1
    if(n == nrow(inputmat)){
      break
    }
    
    srctmp<-rep(i,length(inputmat[n:nrow(inputmat)]))
    dftmp<-inputmat[i,]
    
    
    targettmp<-row.names(as.data.frame(dftmp))[n:nrow(inputmat)]
    
    cortmp<-inputmat[c(n:nrow(inputmat)),i]
    
    
    dataframetmp<-data.frame("Source"=srctmp,"Target"=targettmp,"cor"=cortmp)
    
    dataframetmp<-dataframetmp[dataframetmp$cor> cutoff,]
    
    finaldf<-rbind(finaldf,dataframetmp)
  }
  
  
  finaldf$SourceSym<-mapIds(org.Hs.eg.db,keys=gsub("\\..*","",finaldf$Source),column="SYMBOL",keytype="ENSEMBL",multivals="first")
  
  finaldf$TargetSym<-mapIds(org.Hs.eg.db,keys=gsub("\\..*","",finaldf$Target),column="SYMBOL",keytype="ENSEMBL",multivals="first")
  
  #Check for genes with no Symbol and replace with its ENS identifier 
  for (i in c(1:length(finaldf$SourceSym))){
    
    if (is.na(finaldf$SourceSym[i])){
      
      finaldf$SourceSym[i]<-as.character(gsub("\\..*","",finaldf$Source[i]))
      
    }
    if (is.na(finaldf$TargetSym[i])) {
      finaldf$TargetSym[i]<-as.character(gsub("\\..*","",finaldf$Target[i]))
    }
    
  }
  
  #Check if genes are connecting to itself (b/c ens maps to same SYMBOL) and remove this row
  removerows<-c()
  for(i in c(1:nrow(finaldf))){
    if (finaldf$SourceSym[i] == finaldf$TargetSym[i]){
      removerows<-c(removerows,i)
      next
    }
  }
  if(length(removerows)>0){
  finaldf<-finaldf[-c(removerows),]
  }
  
  #cytodf<-data.frame("Source"=finaldf$Source,"Target"=finaldf$Target,"cor"=finaldf$cor)
  cytodf<-data.frame("Source"=finaldf$SourceSym,"Target"=finaldf$TargetSym,"cor"=finaldf$cor)
  write.csv(cytodf,file=paste("Cyto_networks/",col_mod,"_cyto.txt",sep=""),quote=FALSE,row.names=FALSE)
  
  return(cytodf)
}