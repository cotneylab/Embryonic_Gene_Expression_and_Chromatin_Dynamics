#Expecting input of dataframe of the node files (see table S6 from paper)


WGCNA_random<-function(){
  
  tmpshuffle<-sample(c(1:26122))
  tmpranddf<-wgcna_genespermod
  tmpranddf$color<-wgcna_genespermod$color[tmpshuffle]
  
  return(tmpranddf)
}

permutationfunc<-function(n){
  obs<-c()
  
  for(i in c(1:n)){
   scrambled_wgcna<-WGCNA_random()
   enrichedmod_count<-STRING_db_enrichment(scrambled_wgcna) 
   obs<-c(obs,enrichedmod_count)
   print(paste0("Finished ", i," of ",n,"  permutations."))
  }

  return(obs)
}

STRING_db_enrichment<-function(wgcna_df){
df<-data.frame(NULL)
for(i in unique(wgcna_df$color)){
  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0.4, input_directory="" )
  tmpmod<-subset(wgcna_df,color == i)
  tmpmod_sub<-subset(tmpmod, EHSE == 1 | HUB == 1)
  if(nrow(tmpmod_sub )> 500){
    tmpmod_sub<- tmpmod_sub[1:500,]
  }
  
  tmpmod_mapped<-string_db$map(tmpmod_sub,"ENSID",removeUnmappedRows = TRUE)
  d<-string_db$get_ppi_enrichment(tmpmod_mapped$STRING_id)
  
  r<-data.frame("Enrichment"=d$enrichment,"Lambda"=d$lambda)
  row.names(r)<-i
  
  df<-rbind(df,r)
  
}
df_sub<-subset(df, df$Enrichment <= 0.00172)

return(nrow(df_sub))
}


##To run 
#library(STRINGdb)
#source("stringdb.R")
#load("wgcna_genespermod.Rdata")
#observed<-permutationfunc(20)
#save(observed,file="observed1_ppi_enrich_ct.Rdata")


