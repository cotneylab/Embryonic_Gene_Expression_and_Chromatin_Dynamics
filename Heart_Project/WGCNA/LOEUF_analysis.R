

get_dec_assignment<-function(df){
decli<-c()
for(i in as.character(df$ENSID)){
  
  if( i %in% as.character(LOEUF[[1]]$V2)){
    dec<-"dec1"
  }
  else if( i %in% as.character(LOEUF[[2]]$V2)){
    dec<-"dec2"
  }
  else if( i %in% as.character(LOEUF[[3]]$V2)){
    dec<-"dec3"
  }
  else if( i %in% as.character(LOEUF[[4]]$V2)){
    dec<-"dec4"
  }
  else if( i %in% as.character(LOEUF[[5]]$V2)){
    dec<-"dec5"
  }
  else if( i %in% as.character(LOEUF[[6]]$V2)){
    dec<-"dec6"
  }
  else if( i %in% as.character(LOEUF[[7]]$V2)){
    dec<-"dec7"
  }
  else if( i %in% as.character(LOEUF[[8]]$V2)){
    dec<-"dec8"
  }
  else if( i %in% as.character(LOEUF[[9]]$V2)){
    dec<-"dec9"
  }
  else if( i %in% as.character(LOEUF[[10]]$V2)){
    dec<-"dec10"
  }
  else{
    dec<-"none"
  }
  decli<-c(decli,dec)
  
}
return(decli)
}

perm<-function(df){
  deciles<-c("dec1","dec2","dec3","dec4","dec5","dec6","dec7","dec8","dec9","dec10","none")
  
  
df_rand<-as.data.frame(matrix(ncol=11,nrow=1000))
colnames(df_rand)<-deciles
for(i in c(1:1000)){
  #round(nrow(df)*0.1,0)
  r<-sample(c(1:nrow(df)),2617)
  df_sub<-df[r,]
  
  tmprow<-table(df_sub$dec)

  if(length(tmprow) < 11){
    
    for( j in c(1:length(tmprow))){
      df_rand[i,names(tmprow)[j]]<-unname(tmprow)[j]    
    }
  
    n<-which(!deciles %in% names(tmprow))
    df_rand[i,n]<-0
    
  }
  else{
    df_rand[i,names(tmprow)]<-unname(tmprow)
  }
  
}

return(df_rand)
}


get_LOEUF_permod<-function(mod){
  print(mod)
  deciles<-c("dec1","dec2","dec3","dec4","dec5","dec6","dec7","dec8","dec9","dec10","none")
  currentmod<-read.csv(paste0("/home/tara/R/Human_Heart/WGCNA/Sub_Network/NodeFiles_Rev2/",mod,"_pli_Nodes.csv"),sep=",",header=TRUE)
  
  hub<-subset(currentmod,currentmod$HUB == 1)
  nonhub<-subset(currentmod,currentmod$HUB == 0)
  hub_dec<-get_dec_assignment(hub)
  nonhub_dec<-get_dec_assignment(nonhub)
  hub$dec<-hub_dec
  nonhub$dec<-nonhub_dec
  
  df_nonhub_rand<-perm(nonhub)
 
  med_rand<-c()
  sd_rand<-c()
  for(i in c(1:11)){
    med_rand<-c(med_rand,median(df_nonhub_rand[,i]))
    sd_rand<-c(sd_rand,round(sd(df_nonhub_rand[,i]),0))
  }
  
  hub_counts<-as.data.frame(matrix(ncol=11,nrow=1))
  colnames(hub_counts)<-deciles
  if(length(table(hub$dec)) < 11){
    
    for( i in c(1:length(table(hub$dec)))){
      hub_counts[1,names(table(hub$dec))[i]]<-unname(table(hub$dec))[i]    
    }
    
    n<-which(!deciles %in% names(table(hub$dec)))
    hub_counts[1,n]<-0
    
  }
  else{
    hub_counts[1,names(table(hub$dec))]<-unname(table(hub$dec))
  }

  plothub<-data.frame("Bins"=deciles,"ct"=as.numeric(as.vector(hub_counts[1,])),"sd"=c(rep(NA,11)),"type"=c(paste0(mod,"_Hub")))
  plotnonhub<-data.frame("Bins"=deciles,"ct"=med_rand,"sd"=sd_rand,"type"=c("random_nonHub"))
  
  fullplotdat<-rbind(plothub,plotnonhub)
  fullplotdat$Bins<-as.factor(fullplotdat$Bins)
  fullplotdat$Bins<-factor(fullplotdat$Bins,levels=c("dec1","dec2","dec3","dec4","dec5","dec6","dec7","dec8","dec9","dec10","none"))
  
  p<-ggplot(fullplotdat,aes(x=Bins,y=ct,fill=type))+geom_bar(stat="identity",position="dodge")+scale_fill_manual(values=c(mod,"gray71"))+xlab("")+ylab("Number of genes")+theme(legend.title=element_blank())+geom_errorbar(aes(ymin=ct-sd, ymax=ct+ sd),position=position_dodge(.9),width=0.25)
 return(p) 
}

