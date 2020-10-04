
##Function to get initial modules 
#' @param minMod 
#' @param dpsp
#' @param gTree


dynamicModsfunc<-function(minMod=100,dpsp=2,gTree=geneTree,DT=dissTOM){

	dynamicMods<-cutreeDynamic(dendro=gTree, distM=DT,
			deepSplit=dpsp,pamRespectsDendro=FALSE,
			minClusterSize=minMod)

	dynamicColors=labels2colors(dynamicMods)

	return(dynamicColors)

	
}



##Function to get dendrogram of module eigengenes. 
#' @param Expr the data frame of genes and expression values for samples. Genes are columns.
#' @param dynamicColors the module colors obtained from cutreeDynamic
#' @param deepsplit is the integer from 0 to 4 used in cutreeDynamic
#' @param minmodule is the minModuleSize used in cutreeDynamic

MEtree<-function(Expr,dynamicColors,deepsplit=NULL,minmodule=NULL){
	
	name=paste("module_eigen_dp",deepsplit,"_mm",minmodule,sep="")
	
	#Calculate eigengenes
	MEList<- moduleEigengenes(Expr,colors=dynamicColors)
	MEs<-MEList$eigengenes
	#Calculate dissimilarity of module eigengenes
	MEDiss<-1-cor(MEs)
	#Cluster module eigengenes 
	METree=hclust(as.dist(MEDiss),method='average')
	
	#output the dendrogram 
	pdf(paste("Figures/",name,".pdf",sep=""))
	plot(METree,main=name,xlab="",sub="")
	dev.off()
	return(message("Done"))
	
}


mergemods<-function(Expr,gTree,dynamicColors,MEDissThres,deepsplit,minmodule){
	name=paste("GenedendMerge_dp",deepsplit,"_mm",minmodule,sep="")
	merge=mergeCloseModules(Expr,dynamicColors,cutHeight=MEDissThres,verbose=3)
	mergedColors=merge$colors
	mergedMEs=merge$newMEs
	
	#output the dendrogram with merged modules
	pdf(paste("Figures/",name,".pdf",sep=""))	
	plotDendroAndColors(gTree,cbind(dynamicColors,mergedColors),
	c("Dynamic Tree Cut", "Merged dyanmic"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
	dev.off()
	returnList<-list("mergedCol"=mergedColors,"mergedME"=mergedMEs)
        return(returnList)
        return(message("Done"))
        

}


