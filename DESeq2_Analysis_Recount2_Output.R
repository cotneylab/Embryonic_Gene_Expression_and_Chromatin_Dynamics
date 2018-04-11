## This file gives an example with explanations for how to perform 
# some DESeq2 analysis from Recount2 count tables. 

##Libraries Required 

library(sva)
library(recount)
library(DESeq2)



## Assume the file loads into an object called rse_gene
load("[your recount2 file]")

## Necessary to scale the counts from recount2 object
rse_gene <- scale_counts(rse_gene)

## Create the DESeq object 
dds <- DESeqDataSet(rse_gene, design=~sample)

## Filter
dds<-dds[rowSums(counts(dds))>1,]
## Determine how many genes were filtered out
nrow(dds) 

## Calculate SV variables 
## prep data 
dds<- estimateSizeFactors(dds) 

## Create Model 
mod <- model.matrix( ~sample, colData(dds)) 
mod0 <- model.matrix(~ 1, colData(dds))
 
##Extract just the counts from DESeq object 
dat<- counts(dds,normalized=TRUE)
##Perform svseq 
svseq <- svaseq(dat,mod,mod0)

##Add your SVs to your DESeq object

##If there are not many to add just do it manually: 

ddssva<-dds 
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- ~SV1 + sample

## For when SVs are alot:

for(i in 1:svseq$n.sv) {
	currentSV<- paste("SV", i, sep="")
	ddssva[[currentSV]] <- svseq$sv[,i]
	}
##Then create text of the new formula

form<- "~sample"
formtemp<- c(rep("SV",svseq$n.sv))
for(i in 1:svseq$n.sv){
	formtemp[i] <- paste(formtemp[i],i,sep='')
	form<- paste(form, formtemp[i], sep="+")
	}

## Modify the design to add in your new SVs 
design(ddssva) <- ~sample + SV1 + SV2 + ...SVn

## Run DESeq2 analysis 

ddssva <- DESeq(ddssva)







