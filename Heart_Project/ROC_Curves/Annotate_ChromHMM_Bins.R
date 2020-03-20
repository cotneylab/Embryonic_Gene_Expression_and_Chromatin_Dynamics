#The posterior probabilities from ChromHMM are in 200bp bins that are not annotated.
#This function adds the coordinates. Need to run for every chromosome. Used Hg19 sizes.

sep_chr_bins<-function(){
#Just storing where I did this. Files for all 23 chr in folder already.
n="chr21" #example!!
file<-grep(x=list.files("CS13"),pattern=paste0(n,"_posterior.txt"),value=TRUE)
chr_tmp<-read.table(paste0("CS13/",file),skip=1,sep="\t",header=TRUE)

samp_from<-seq(from=1, to=48129601, by=200)
samp_to<-seq(from=200, to=48129800, by=200)

#check to make sure lengths match the chromHMM file
isTRUE(length(samp_from)==length(samp_to) &&   
  length(samp_from)==nrow(chr_tmp)) 

chr<-c(rep(n,length(samp_from)))

bin_bed<-data.frame("chr"=chr,"From"=samp_from,"To"=samp_to)

save(bin_bed,file=paste0(n,"_binbed.Rdata"))
}