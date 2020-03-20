library(GenomicRanges)

TN<-read.table("~/ANALYSIS/ChIP-seq/Human_Heart/ROC/Vista_positive_no_heart_hg19.bed")
TP<-read.table("~/ANALYSIS/ChIP-seq/Human_Heart/ROC/Vista_positive_heart_hg19.bed")

colnames(TP)<-c("chr","start","end")
colnames(TN)<-c("chr","start","end")

TP_gr<-makeGRangesFromDataFrame(TP)

TN_gr<-makeGRangesFromDataFrame(TN)

peaks_gr<-makeGRangesFromDataFrame(peaks)

TP_overlap<-findOverlaps(TP_gr,peaks_gr)

TN_overlap<-findOverlaps(TN_gr,peaks_gr)

TP_overlap<-as.data.frame(TP_overlap)

TN_overlap<-as.data.frame(TN_overlap)

source("ROC_calc.R") (in ROC_emerge directory but see below aswell)

TN_ROC<-pred_per_peak(TN,TN_overlap)

TP_ROC<-pred_per_peak(TP,TP_overlap)

emerge_ROC<-data.frame("label"=c(rep(0,846),rep(1,281)),"prediction"=c(TN_ROC,TP_ROC))

save(emerge_ROC,file="emerge_ROC_data.Rdata")

#This script will sum the prediction value from the emerge peak file
#for each TP or TN segment

pred_per_peak<-function(test_segments,overlaps){
#test_segments: bed file of True positives or True negatives
#overlaps: dataframe version of output from overlaps function of
#GenomicRanges

pred_per_peak<-c()
for(i in c(1:nrow(test_segments))){

tmp<-subset(overlaps, overlaps$queryHits == i)
if(nrow(tmp) == 0){
tmp_pred = 0
pred_per_peak<-c(pred_per_peak,tmp_pred)
next
}

tmp_pred<-sum(peaks$pred[tmp$subjectHits])
#peaks is a bed file where fourth column is the score

pred_per_peak<-c(pred_per_peak,tmp_pred)
}
return(pred_per_peak)
}

