
#function called in get_avg_pertiss.sh
#This function will get the avg expression per tissue for each of the genes from the GREAT analysis 
#Expecting a dataframe with the ensembl IDs for the inquired genes
#Expecting an empty dataframe to fill up
#Expecting to beable to pull from recount2 GTEx norm counts matrix
#Returns dataframe with average expression values

avg_pertiss<-function(avgXpress){
    
for(n in colnames(avgXpress)){  
print(n)
filename<-paste0("~/ANALYSIS/RNA-seq/Human_Heart/GTEx_compare/GTEX_norm_cts/","G
tex_",n,"norm_cts.Rdata")
load(filename)
tmp<-as.data.frame(tmp)

df<-NULL
for(i in c(1:nrow(enh_cts_sym_ens))){
  
  indx<-which(enh_cts_sym_ens$ENS[i] == gsub("\\..*","",row.names(tmp)))
  
  df<-rbind(df,tmp[indx,])
  
}
avgpergene<-rowMeans(df)

avgXpress[,n]<-avgpergene

}

return(avgXpress)
}

