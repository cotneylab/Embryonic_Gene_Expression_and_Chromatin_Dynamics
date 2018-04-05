#Function to compute z-score of a row
#Input: counts ~ dataframe of raw values
#Output: zscore ~dataframe of corresponding z-score


zscore <- function(counts) {
	rownum <- dim(counts)[1]
	colnum <- dim(counts)[2]
	zscore <- data.frame(matrix(NA,nrow=rownum,ncol=colnum))
	rownames(zscore) <- rownames(counts)
	colnames(zscore)<-colnames(counts) 

	for(i in 1:dim(counts)[1]) {
		zscoretemp <- ( counts[i,] - mean(counts[i,]) ) / sd(counts[i,])    
		zscore[i,]<-zscoretemp	

	}

	return(zscore)
}
