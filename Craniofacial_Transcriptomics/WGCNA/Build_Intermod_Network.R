library(ggraph)
library(igraph)

get_colors<-function(dat){
  #Get the module colors
  colvec<-c()
  
  for(i in row.names(dat)){  
    
    tempcol<-substr(i, start=3, stop=nchar(i))
    
    colvec<-c(colvec,tempcol)
    
  }
  return(colvec)
  
}

get_MEcol<-function(dat){
  
  n<-c()
  for (i in row.names(dat)){
    tmp<-paste("ME",i,sep="")
    n<-c(n,tmp)
  }

  return(n)
}


makeSegments<-function(cor_ME,mod_points,segmentdf,modules){
  
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
  
  colnames(endpts)<-c("V1e","V2e")
  endpts$mod<-row.names(endpts)
  tmpli<-mod_points[currentmod,]
  
  startpts<-as.data.frame(lapply(tmpli,rep,(nrow(cor_ME_tmp)-1)))

  segments<-cbind(startpts,endpts)
 
  #segments$to<-c(row.names(startpts),row.names(endpts))
  
  segments$cor<-cor_ME_tmp[,currentmod][2:nrow(cor_ME_tmp)]
  
  segments$node<-c(rep(currentmod,nrow(segments)))
  
  segmentdf<-rbind(segmentdf,segments)
  
  return(makeSegments(cor_ME=cor_ME,mod_points=mod_points,segmentdf=segmentdf,
                      modules=modules[2:length(modules)]))
}



