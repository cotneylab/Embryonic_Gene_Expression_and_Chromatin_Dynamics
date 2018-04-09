#Function to extract GO term or description from the "Term" column
#you get from david and collect in vector.
#Make sure your input column is converted to  character vector
#as.character(<termcolumn>)
#n = 1 for when you just want GO number, n=2 for description.

GOcolumn<-function(termcolumn,n=1){
  
  splitvector<-strsplit(termcolumn,split="~")  
  GOcol<-character()
  
  for(i in 1:length(splitvector)) {
  
  
    tempGO<- splitvector[[i]][n] 
    GOcol<-c(GOcol,tempGO)
  
  
  }
  return(GOcol)
}

