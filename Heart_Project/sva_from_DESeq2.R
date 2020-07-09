##General code to get SVA corrected counts from DESeq2 objects

library(DESeq2)

##load your ddssva and SVs
load("<ddssva filename>")
load("<svseq filename>")

##Extract the counts
dat<-counts(ddssva,normalized=TRUE)

## turn counts dataframe into type matrix
dat<-as.matrix(dat)

##transpose of count matrix
Y<- t(dat)

## Extract just the SVs
W<- svseq$sv

##complicated matrix math
alpha <- solve(t(W) %*% W) %*% t(W) %*% Y

##Add sva batch corrected count values to the svseq object
svseq$corrected<- t(Y- W %*% alpha)

##save your svseq object
save(svseq,file="HumanHeart_rnaseq_SVcorrect.Rdata")
