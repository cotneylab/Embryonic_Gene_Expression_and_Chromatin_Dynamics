#This script will sum the prediction value from the emerge bedgraph file for each peak  ~84443 peaks total 

overlaps<-findOverlaps(peaks_gr,bedgraph_gr)
#xx_gr is genomicRanges object of bed file for peaks or bedgraph

save(overlaps,file="emerge_peak_bedgraph_overlap.Rdata")

library(GenomicRanges)
load("emerge_pred_track.Rdata")

load("emerge_peak_bedgraph_overlap.Rdata") 
overlaps<-as.data.frame(overlaps)

pred_per_peak<-c()
for(i in c(1:84443)){

tmp<-subset(overlaps, overlaps$queryHits == i)

tmp_pred<-sum(bedgraph$V4[tmp$subjectHits])

pred_per_peak<-c(pred_per_peak,tmp_pred)
}

save(pred_per_peak,file="predperpeak.Rdata")
