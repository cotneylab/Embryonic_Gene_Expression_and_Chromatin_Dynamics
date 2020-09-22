#Tiss is a list of tissues
#ginigenes are the dataframe of the gini-genes with tissue assignment
#li is the list of genes you are looking for enrichment


permutationfunc<-function(tiss){
  obs<-NULL
  
  for(i in c(1:1000)){
    tmpshuffle<-sample(c(1:nrow(ginigenes)))
    tmpranddf<-ginigenes
    tmpranddf$Tissue<-ginigenes$Tissue[tmpshuffle]
    obs_r<-NULL
    for(j in tiss){
      tissgenes<-subset(tmpranddf,Tissue == j)
      
      overlap<-length(intersect(tissgenes$ENS,li))
      obs_r<-c(obs_r,overlap)
    }
   
    obs<-rbind(obs,obs_r)
  }
  return(obs)
}

#df is dataframe

get_sd<-function(df){
  sd<-c()
for(i in c(1:34)){
  sd_tmp<-round(sd(df[,i]),0)
  sd<-c(sd,sd_tmp)
}
  return(sd)
}

