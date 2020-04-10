##Be sure to have logcts (see how to setup)

#load("HumanHeart_rnaseq_dds_rlog.Rdata")
#load("Human_Heart_WGCNA_mods.Rdata")
#logcts<-assay(dds_filt_log)
#logcts<-as.data.frame(logcts)
#logcts$col<-merged_mods$mergedCol
#row.names(logcts)<-gsub("\\..*","",row.names(logcts))
#modcols<-as.list(unique(logcts$col))


##Input a list of colors for "mod"
##Input a character vector of genes you want to check enrichment for as "li"

Enrichment_count<-function(mod){
  li<-as.character(Enh_cts_sym_ENS_filt$ENS)

  modgenes<-subset(logcts, logcts$col == mod)
  
  overlap<-length(intersect(row.names(modgenes),li))
  
  return(overlap)
}

Enrichment_count_random<-function(mod){
  li<-repbound_ens2
  
  tmpshuffle<-sample(c(1:26122))
  tmpranddf<-logcts
  tmpranddf$col<-logcts$col[tmpshuffle]
  
  modgenes<-subset(tmpranddf, tmpranddf$col == mod)
  
  overlap<-length(intersect(row.names(modgenes),li))
  
  return(overlap)
}

permutationfunc<-function(){
  obs<-NULL
  
  for(i in c(1:1000)){
    random_Enrichment<-sapply(modcols,Enrichment_count_random,simplify=FALSE)
    obs<-rbind(obs,random_Enrichment)
  }
  return(obs)
}

pval_func<-function(actual,observed){
  
  #last step
  
  pval<-c()
  
  oddsratio<-c()
  
  
  for(i in c(1:length(actual))){
    
    tmppval<-(sum(observed[,i] > actual[[i]])+1) / (1000+1)
    
    tmpor<- actual[[i]] / median(unlist(observed[,i]))
    
    pval<-c(pval,tmppval)
    
    oddsratio<-c(oddsratio,tmpor)
    
  }
  df<-data.frame("mod"=unique(logcts$col),"pval"=pval,"OR"=oddsratio)
  return(df)
  
}

##Example of how to use these functions
##NOTE: Need to change the gene list in the Enrichment functions

## To get the actual values

#actual<-sapply(modcols,Enrichment_count,simplify  = FALSE)

##To get the random values (1000) which calls enrichment_count_random

#observed<-permutationfunc()

##Get the pvals and OR 

#HSE_enrichment<-pval_func(actual,observed)

#HSE_enrichment<-HSE_enrichment[order(HSE_enrichment$pval),]





