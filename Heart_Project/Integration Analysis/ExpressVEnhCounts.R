>Scatterplot of gene expression vs Embryonic Heart Specific Enhancers 
>Extract the normalized counts per tissues 
#Liver, Spleen, Thyroid, Blood, Brain

load("Gtex_emheart_DESeq2_norm_cts.Rdata")

m<-read.csv("Recount2_Gtex_Heart_metadata_map.csv",header=TRUE)

brainindx<-which(m$smts == "Brain")
brain<-dat[,brainindx]

liverindx<- which(m$smts == "Liver")
liver<-dat[,liverindx]

spleenindx<-which(m$smts == "Spleen")
spleen<-dat[,spleenindx]

thyroidindx<- which(m$smts == "Thyroid")
thyroid<-dat[,thyroidindx]

bloodindx<-which(m$smts == "Blood")
blood<-dat[,bloodindx]

save(blood,file="Gtex_blood_norm_cts.Rdata")
save(spleen,file="Gtex_spleen_norm_cts.Rdata")
save(liver,file="Gtex_liver_norm_cts.Rdata")
save(thyroid,file="Gtex_thyroid_norm_cts.Rdata")
save(brain,file="Gtex_brain_norm_cts.Rdata")
save(emheart,file="Gtex_emheart_norm_cts.Rdata")


>Get list of EHSE genes and their enh counts 
#Genes and enhancers from Jen Note contains 8208 genes 
jen_genes<-read.csv("great_enhancers_assigned_to_gene.tsv",sep="\t")
#List of ENS to SYM from great 
m<-read.csv("hg19.great3.0.genes.txt",sep="\t",header=FALSE)
#Build Dataframe #Note, 12 genes from jens list not in the great map (all 43xxx genes). Total matrix is 8196 genes 
df<-NULL
for(i in c(1:nrow(jen_genes))){
  
  indx<-which(jen_genes$X..GREAT.version.3.0.0[i] == as.character(m$V5))
  if(length(indx) > 1){
    indx <- indx[1]
  }
  if(length(indx) == 0){
    print(jen_genes$X..GREAT.version.3.0.0[i])
    next
  }
  
  tempdf<-cbind(jen_genes[i,], data.frame("ENS"=as.character(m$V1)[indx]))
  
  df<-rbind(df,tempdf)
}

#Get the counts of enhancers for each gene 
enhancer_counts<-c()
for(i in df$Species.assembly..hg19){
  
  temp<-unlist(strsplit(as.character(i),split=",",fixed=TRUE))
  
  enhancer_counts<-c(enhancer_counts,length(temp))
  
}
colnames(df)<-c("SYM","Enhloc","ENS","Enhcts")
save(Enh_cts_sym_ENS,file="Great_Enh_cts_ENS_Sym.Rdata")




load("Great_Enh_cts_ENS_Sym.Rdata")
#Get the rows that match in Recount2 genes (total of 8113)
df<-NULL
for(i in c(1:nrow(Enh_cts_sym_ENS))){
  
  indx<-which(Enh_cts_sym_ENS$ENS[i] == gsub("\\..*","",row.names(emheart)))
  
  df<-rbind(df,emheart[indx,])
  
}

#remove the PAR_Y genes (now have a total of 8111 genes) 
df<-df[-c(519,7830),]
#Create subset for the all of the tissues 
EHS_spleen<-as.data.frame(spleen[row.names(df),])
EHS_liver<-as.data.frame(liver[row.names(df),])
EHS_brain<-as.data.frame(brain[row.names(df),])
EHS_blood<-as.data.frame(blood[row.names(df),])
EHS_emheart<-as.data.frame(emheart[row.names(df),])
EHS_thyroid<-as.data.frame(thyroid[row.names(df),])


#Create subset of the zero enhancer target genes 
EHS_emheart_zo<-as.data.frame(subset(emheart, !(row.names(emheart) %in% row.names(df))))
EHS_brain_zo<-as.data.frame(subset(brain, !(row.names(brain) %in% row.names(df))))
EHS_blood_zo<-as.data.frame(subset(blood, !(row.names(blood) %in% row.names(df))))
EHS_liver_zo<-as.data.frame(subset(liver, !(row.names(liver) %in% row.names(df))))
EHS_spleen_zo<-as.data.frame(subset(spleen, !(row.names(spleen) %in% row.names(df))))
EHS_thyroid_zo<-as.data.frame(subset(thyroid, !(row.names(thyroid) %in% row.names(df))))




#filter the zeros for low counts
EHS_zero_filt<-cbind(EHS_emheart_zo,EHS_blood_zo,EHS_brain_zo,EHS_liver_zo,EHS_spleen_zo,EHS_thyroid_zo)
EHS_zero_filt<-EHS_zero_filt[rowSums(EHS_zero_filt)>100,]
EHS_emheart_zo<-EHS_emheart_zo[row.names(EHS_zero_filt),]
EHS_blood_zo<-EHS_blood_zo[row.names(EHS_zero_filt),]
EHS_brain_zo<-EHS_brain_zo[row.names(EHS_zero_filt),]
EHS_liver_zo<-EHS_liver_zo[row.names(EHS_zero_filt),]
EHS_spleen_zo<-EHS_spleen_zo[row.names(EHS_zero_filt),]
EHS_thyroid_zo<-EHS_thyroid_zo[row.names(EHS_zero_filt),]

#Setup plots 
plotdf_zo<-data.frame("AvgX"=c(rowMeans(EHS_emheart_zo),rowMeans(EHS_blood_zo),rowMeans(EHS_brain_zo),rowMeans(EHS_liver_zo),rowMeans(EHS_spleen_zo),rowMeans(EHS_thyroid_zo)),
"Tissue"=c(rep("Emheart",43425),rep("Blood",43425),rep("Brain",43425),rep("Liver",43425),rep("Spleen",43425),rep("Thyroid",43425)))
save(plotdf_zo,file="plotting_zero_enh.Rdata")
#Get log2 values 
plotdf_zo$log2AvgX<-log(plotdf_zo$AvgX+1,2)
plotdf_zo$ENS<-c(rep("dummy",260550))
plotdf_zo<-plotdf_zo[,c("AvgX","Tissue","ENS","Enhcts","log2AvgX")]


plotdf<-data.frame("AvgX"=c(rowMeans(EHS_emheart),rowMeans(EHS_blood),rowMeans(EHS_brain),rowMeans(EHS_liver),rowMeans(EHS_spleen),rowMeans(EHS_thyroid)),
"Tissue"=c(rep("Emheart",8111),rep("Blood",8111),rep("Brain",8111),rep("Liver",8111),rep("Spleen",8111),rep("Thyroid",8111)))
save(plotdf,file="plotting_Enhtargets.Rdata")

#add gene names 
plotdf$ENS<-gsub("\\..*","",row.names(df))
#Add enh counts 
Enh_cts_sym_ENS_filt<-Enh_cts_sym_ENS[gsub("\\..*","",row.names(df)),]
plotdf$Enhcts<-Enh_cts_sym_ENS_filt$Enhcts
#Get log2 values 
plotdf$log2AvgX<-log(plotdf$AvgX+1,2)

save(plotdf_zo,file="plotting_zero_enh_final.Rdata")
save(plotdf,file="plotting_Enhtargets_final.Rdata")
plotdf_final<-rbind(plotdf,plotdf_zo)

#Making Boxplots 
boxpl_emheart<-subset(plotdf_final, Tissue =="Emheart")
boxpl_blood<-subset(plotdf_final, Tissue =="Blood")

boxpl<-rbind(boxpl_emheart,boxplot_blood)

vec<-c()
for(i in as.numeric(boxpl$Enhcts)){
  
  if(i == 0){
    temp<-"zero"
  }
  else if(i >=1 && i <25){
    temp<-"oneto25"
  }
  else if(i > 25) {
    temp<-"over25"
  }
  vec<-c(vec,temp)
}

boxpl$bin<-vec
#CRITICAL PIECE OF CODE TO MAKE LIFE EASIER FOR PLOTTING#
boxpl$Tissue<-factor(boxpl$Tissue,levels=levels(boxpl$Tissue)[c(3,1,2,4,5,6)])

 pbx<-ggplot(data=boxpl,aes(x=bin, y=log2AvgX, fill=Tissue))+geom_violin(draw_quantiles = c(0.25,0.5,0.75))+stat_compare_means(aes(group=Tissue),label="p.format")+scale_fill_manual(values=c("orangered2","gray71"))

#Make the significance plots:

ggplot(data=second,aes(x=Tissue,y=log2avg,fill=Tissue))+geom_violin(draw_quantiles = c(0.25,0.5,0.75))+scale_fill_manual(values=c("orangered2",rep("grey71",31)))+stat_compare_means(ref.group = c("Embryonic Heart"),label="p.signif")+theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1))


