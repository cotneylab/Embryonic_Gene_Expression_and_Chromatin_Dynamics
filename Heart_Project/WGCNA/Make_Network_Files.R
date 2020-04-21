
library(AnnotationDbi)
library(org.Hs.eg.db)
load("gene2gini.Rdata")
load("Enh_cts_sym_ens.Rdata")
load("known_chd.Rdata")
load("chd_rare_exome.Rdata")
load("nkx25_act_bound_ens.Rdata")

load("~/R/Human_Heart/WGCNA/Human_Heart_WGCNA_mods.Rdata")
load("~/R/Human_Heart/WGCNA/HumanHeart_rnaseq_dds_rlog.Rdata")
load("~/R/Human_Heart/WGCNA/Intramod_Connectivity.Rdata")

logcts<-assay(dds_filt_log)
logcts<-as.data.frame(logcts)
logcts$col<-merged_mods$mergedCol
row.names(logcts)<-gsub("\\..*","",row.names(logcts))

modcols<-as.list(unique(logcts$col))

emheartginipt5<-subset(gene2gininew, Max_Tissue == "Embryonic Heart")
emheartginipt5<-subset(emheartginipt5, Gini >= 0.5)


ENS2SYM<-function(vec){
  ens2sym<-mapIds(org.Hs.eg.db,keys=gsub("\\..*","",vec),column="SYMBOL",keytype="ENSEMBL",multivals="first")
 
   #For those without a symbol, just keep ENS name 
  for(i in c(1:length(ens2sym))){
    if(is.na(ens2sym[i])){
      ens2sym[i]<-as.character(vec[i])
    }
  }
  return(unname(ens2sym))
}

make_attribute_list<-function(nodeID,genelist){
  attribute_value<-c()
  attribute_name<-c()
  
  for(i in nodeID){
  
  if(i %in% genelist){
    
    temp_attribute_value<-c(1)
    temp_attribute_name<-c(i)
  }
  else{
    temp_attribute_value<-c(0)
    temp_attribute_name<-c("NA")
  }
  attribute_value<-c(attribute_value,temp_attribute_value)
  attribute_name<-c(attribute_name,temp_attribute_name)

  }
  
  return(list("val"=attribute_value,"name"=attribute_name))

}

node_files<-function(module){
  
modgenes<-subset(logcts, logcts$col == module)
ID<-row.names(modgenes)

#Get top 10% intra-connected genes to be used as "Hub" genes 
Intracon<-Alldegrees1[row.names(modgenes),]
Intracon<-Intracon[order(-Intracon$kWithin),]
Hub<-gsub("\\..*","",row.names(Intracon)[1:round(nrow(Intracon)*.1,0)])

#Convert all ENS to gene symbols 
Node_sym<-ENS2SYM(ID)
Gini_gene_sym<-ENS2SYM(emheartginipt5$ENS)
HSE_gene_sym<-ENS2SYM(enh_cts_sym_ens$ENS)
Hub_sym<-ENS2SYM(Hub)
Rare_Disease_sym<-ENS2SYM(rare_chd$ens)
Known_Disease_sym<-ENS2SYM(known_chd$ENS)
nkx25bound_sym<-ENS2SYM(nkx25Actbound_ens)

gini<-make_attribute_list(Node_sym,Gini_gene_sym)
hse<-make_attribute_list(Node_sym,HSE_gene_sym)
hub<-make_attribute_list(Node_sym,Hub_sym)
rare<-make_attribute_list(Node_sym,Rare_Disease_sym)
known<-make_attribute_list(Node_sym,Known_Disease_sym)
nkx25<-make_attribute_list(Node_sym,nkx25bound_sym)

Node_df<-data.frame("ENSID"=ID,"SYMID"=Node_sym, "GINI" = gini$val,"EHSE"=hse$val,"HUB"=hub$val,"Known_Disease"=known$val,
                    "Rare_Disease"=rare$val,"nkx25_bound_act"=nkx25$val,
                    "GINIname"=gini$name,"EHSEname"=hse$name,"HUBname"=hub$name,"Known_Diseasename"=known$name,
                    "Rare_Diseasename"=rare$name,"nkx25_bound_actname"=nkx25$name)


s<-rowSums(Node_df[,3:7])
Node_df$score<-s
Node_df<-Node_df[order(-Node_df$score),]
write.table(Node_df,file=paste("NodeFiles_Rev2/",module,"_Nodes.csv",sep=""),sep=",",quote = FALSE,row.names=FALSE)
print(paste0("finished ",module))
}

############################################
#To run:

#mcols<-unique(merged_mods$mergedCol)
#for(i in mcols){
  
#  node_files(i)
  
#}


