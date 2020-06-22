
#Input df from file MSH2_Windows_Rev2.tsv

SGE_oligos<-function(msh2){

msh2_snvs<-NULL
tmpdf<-NULL
for(j in c(1:nrow(msh2))){

tmpdf<-paste0(as.character(msh2$Pool_F[j]),as.character(msh2$Intron_L[j]),as.character(msh2$Exon[j]),as.character(msh2$Intron_R[j]),as.character(msh2$Pool_R[j]))
bp<-c("A","G","C","T")
nsnv<-as.character(msh2$SNV_pos[j])

pool1<- as.character(msh2$Pool_F[j])
pool2<- as.character(msh2$Pool_R[j])
intron1<- as.character(msh2$Intron_L[j])
intron2<- as.character(msh2$Intron_R[j])
exon<-as.character(msh2$Exon[j])

for(i in c(1:nchar(exon))){
  
  syn_v<-as.numeric(strsplit(as.character(msh2$SNV_pos[j]),",")[[1]])
  
  if( length(syn_v) > 1){
     if(i %in% syn_v){
       next
     }
  }
  else if(i %in%  syn_v){
    next 
  }
  
  snvs<-bp[!grepl(substr(exon,i,i),bp)]
  
  tempexon<-exon
  
  substr(tempexon,i,i)<-snvs[1]
  tempsnv1<-paste(pool1,intron1,tempexon,intron2,pool2,sep="")
  tempexon<-exon
  
  substr(tempexon,i,i)<-snvs[2]
  tempsnv2<-paste(pool1,intron1,tempexon,intron2,pool2,sep="")
  tempexon<-exon
  
  substr(tempexon,i,i)<-snvs[3]
  tempsnv3<-paste(pool1,intron1,tempexon,intron2,pool2,sep="")
  
  tmpdf<-rbind(tmpdf,tempsnv1,tempsnv2,tempsnv3)
  
}
msh2_snvs<-rbind(msh2_snvs,tmpdf)

}

write.table(msh2_snvs,file="MSH2_Exon11_Oligos.tsv",sep="\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
return(msh2_snvs)

}
get_length<-function(){
ct<-c()
for(i in c(1:nrow(windows))){
  
  tmpct<-nchar(as.character(windows[i,5]))+nchar(as.character(windows[i,6]))+nchar(as.character(windows[i,7]))+nchar(as.character(windows[i,8]))+
    nchar(as.character(windows[i,9]))
  
  ct<-c(ct,tmpct)
}
}