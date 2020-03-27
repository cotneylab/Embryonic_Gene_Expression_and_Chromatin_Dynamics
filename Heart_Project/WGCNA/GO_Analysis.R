
parse_mods<-function(dat){
for(mod in unique(merged_mods$mergedCol)){
  
  name<-paste("mod",mod,".RData",sep="")
  modtmp<-dat[which(dat$Module == mod),]
  save(modtmp,file=name)
  
}
}

addlist_david<-function(name){
  print(name)
  
  load(paste("mod",name,".RData",sep=""))
  
  addList(david,gsub("\\..*","",row.names(modtmp)),idType="ENSEMBL_GENE_ID",
          listName=name,listType="Gene")
  
  
}

david_outputs<-function(module_list){
  
  cl<-c(1:length(module_list))
  
  setAnnotationCategories(david, 
                          c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
  for(i in cl){
    
    setCurrentGeneListPosition(david,i)
    setCurrentBackgroundPosition(david,1)
    
    FuncAnnotChart<-getFunctionalAnnotationChart(david)
    getFunctionalAnnotationChartFile(david,file=paste(module_list[i],"FuncAnnotChart.tsv",sep=""))
    
    setCurrentGeneListPosition(david,i)
    setCurrentBackgroundPosition(david,1)
    
    FuncAnnotClust<-getClusterReport(david)
    getClusterReportFile(david,file=paste(module_list[i],"FuncAnnotClust.tsv",sep=""))
  }
  return("Done")
}


make_dotplots<-function(module){

name=paste(module,"FuncAnnotChart.tsv",sep="")

modtmp<-read.csv(name,sep="\t",header=TRUE)

terms<-character()

for (i in modtmp$Term) {
  
  temp<-unlist(strsplit(as.character(i),split="~",fixed=TRUE))[2]
  
  terms<-c(terms,temp)
  
}

modtmp<-modtmp[c("Term","X.","Benjamini")]
modtmp$Term<-terms
modtmp<-modtmp[modtmp$Benjamini<=0.05,]
modtmp<-modtmp[order(modtmp$Benjamini),]
modtmp$Bonferroni<--log(modtmp$Benjamini,2)
modtmp$Term<-factor(modtmp$Term,levels=modtmp$Term[nrow(modtmp):1])
colnames(modtmp)[2]<-c("GenesAssociated")

if(nrow(modtmp) > 0) {

  if( nrow(modtmp) < 25 ) {
pdf(paste(module,"dotplot.pdf"))
p<-ggplot(modtmp,aes(x=Bonferroni,y=Term))
p2<-p+geom_dotplot(binaxis='y',binpositions='all',dotsize=0,binwidth=1.5)+
  geom_point(aes(fill=GenesAssociated),size=5,colour="black",pch=21)+
  scale_fill_gradient2(low="black",mid="gray",high=module)+
  theme(axis.text.y=element_text(size=10)) +ylab(module)
print(p2)
dev.off()
return(paste(module,"module has less than 25 significant terms",sep=" "))
  }
  else if(nrow(modtmp) > 25 ){
  pdf(paste(module,"dotplot.pdf"))
  p<-ggplot(modtmp[1:25,],aes(x=Bonferroni,y=Term))
  p2<-p+geom_dotplot(binaxis='y',binpositions='all',dotsize=0,binwidth=1.5)+
    geom_point(aes(fill=GenesAssociated),size=5,colour="black",pch=21)+
    scale_fill_gradient2(low="black",mid="gray",high=module)+
    theme(axis.text.y=element_text(size=10)) +ylab(module)
  print(p2)
  dev.off()
  return(paste(module,"module has more than 25 significant terms",sep=" "))
  }
}

else{
  return(paste(module,"module has no significant terms",sep=" "))
}

}





