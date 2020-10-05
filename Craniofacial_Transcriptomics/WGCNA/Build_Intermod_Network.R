makeSegments<-function(cor_ME,mod_points,segmentdf=data.frame(NULL),modules){
  
  if(length(modules) == 1){
   
    return(segmentdf)
  }
  
  
  currentmod<-modules[1]
  
  cor_ME_tmp<-subset(cor_ME, cor_ME[,currentmod] >= .5 )
  
  #Goto next module if none meet correlation criteria
  if(nrow(cor_ME_tmp) == 1){
    return(makeSegments(cor_ME=cor_ME,mod_points=mod_points,segmentdf=segmentdf,
                        modules=modules[2:length(modules)]))
  }
  
  cor_ME_tmp<-cor_ME_tmp[order(-cor_ME_tmp[,currentmod]),]
  
  endpts<-mod_points[row.names(cor_ME_tmp)[2:nrow(cor_ME_tmp)],]
  if(nrow(cor_ME_tmp) == 2){
    endpts<-data.frame("V1e"=endpts[1], "V2e"=endpts[2])
  }
  print(endpts)
  
  colnames(endpts)<-c("V1e","V2e")
  endpts<-as.data.frame(endpts)
  endpts$mod<-row.names(endpts)
  tmpli<-mod_points[currentmod,]
  
  startpts<-data.frame("s1"=rep(tmpli[1],nrow(cor_ME_tmp)-1),"s2"=tmpli[2])
  
  segments<-cbind(startpts,endpts)
  
  segments$cor<-cor_ME_tmp[,currentmod][2:nrow(cor_ME_tmp)]
  
  segments$node<-c(rep(currentmod,nrow(segments)))
  
  segmentdf<-rbind(segmentdf,segments)
  
  return(makeSegments(cor_ME=cor_ME,mod_points=mod_points,segmentdf=segmentdf,
                      modules=modules[2:length(modules)]))
}
