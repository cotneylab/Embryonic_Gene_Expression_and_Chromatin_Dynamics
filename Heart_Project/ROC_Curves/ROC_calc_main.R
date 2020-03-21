library("GenomicRanges")
library("dplyr")

source("~/ANALYSIS/ChIP-seq/Human_Heart/ROC/ROC_calc.R")

TP<-read.table("~/ANALYSIS/ChIP-seq/Human_Heart/ROC/Vista_positive_heart_hg19.bed")
TN<-read.table("~/ANALYSIS/ChIP-seq/Human_Heart/ROC/Vista_positive_no_heart_hg19.bed")
colnames(TP)<-c("chr","start","end")
colnames(TN)<-c("chr","start","end")

TP_gr<-makeGRangesFromDataFrame(TP)

TN_gr<-makeGRangesFromDataFrame(TN)


samples_No<-c("12383","12690","12408","14135","12997","14315","12291","12331","12456","12059","11914","12135","12451","12448","11849","12093","12151","12058")
samples_CS<-c(rep("CS13",2),rep("CS14",2),rep("CS16",2),rep("CS17",2),rep("CS18",2),rep("CS19",2),rep("CS20",2),rep("CS21",2),rep("CS23",2))
ROC_data<-vector("list",length(samples_No))

#Loop through for each sample 
for(i in c(1:length(samples_No))){

TP_pred_tmp<-get_overlaps(TP_gr,samples_No[i])

TN_pred_tmp<-get_overlaps(TN_gr,samples_No[i])

df_tmp<-data.frame("labels"=c(rep(0,846),rep(1,281)),"predictions"=c(TN_pred_tmp,TP_pred_tmp))

ROC_data[[i]]<-df_tmp

}

header<-c()
for(i in c(1:length(samples_No))){
header<-c(header,paste0(samples_CS[i],"_",samples_No[i]))
}

names(ROC_data)<-header
save(ROC_data,file="Human_Heart_ROC.Rdata")


