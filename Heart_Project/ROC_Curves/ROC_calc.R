
get_score<-function(overlaps_df,post_prob){

#overlaps: dataframe version of output from overlaps function of GenomicRanges
#post_prob: dataframe containing post prob of each Chrom HMM state.

scores_per_chr<-c()

for(i in unique(overlaps_df$queryHits)){

tmp<-subset(overlaps_df, overlaps_df$queryHits == i)

tot_window<-post_prob[tmp$subjectHits,]

df_sub<-tot_window[,c("E13","E14","E15","E18")]

scores_per_chr<-c(scores_per_chr,mean(rowSums(df_sub)))

}
return(scores_per_chr)
}

get_overlaps<-function(test_segments,sampleNo){
   
  loc<-"/home/FCAM/tyankee/ANALYSIS/ChIP-seq/Human_Heart/ROC/"
  #initialize variables 
  scores<-c()
chrli<-paste0("chr",c(as.character(c(1:22)),"X"))
overlaps<-c()

  for(i in chrli){
  
      load(paste0(loc,grep(x=list.files(loc),pattern=paste0(i,"_binbed.Rdata"),value=TRUE)))
      prob_file<-grep(x=list.files(),pattern=paste0(sampleNo,"_25_",i,"_"),value=TRUE)
      prob_df<-read.table(prob_file,skip=1,sep="\t",header=TRUE)
      prob_df<-cbind(bin_bed,prob_df)
      colnames(prob_df)[1:3]<-c("chr","start","end")
      prob_df_gr<-makeGRangesFromDataFrame(prob_df)
  
    overlap_subset<-findOverlaps(test_segments,prob_df_gr)   
    scores_tmp<-get_score(as.data.frame(overlap_subset),prob_df)

   scores<-c(scores,scores_tmp) 
   
  }

  return(scores)
}
