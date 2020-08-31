library(ineq)

#input: list of tissues
#tissuemap: dataframe of the Recount2 metadata
tissue_colval<-function(input){
  
  return(which(input == as.character(tissuemap$smts)))
  
}


#input:Lists of colnumbers for each sample and tissue
#cts:Dataframe of scaled counts 

Tissue_Average<-function(input){
  tissavg<-c()
  
  for( i in c(1:nrow(cts))){
     tempmean<-round(mean(cts[i,input]),2)
   

    tissavg<-c(tissavg,tempmean)
   
}
  

 return(tissavg)
  
}

#inputdf avg expression matrix, rows are genes and cols are tissue

Tissue_Max<-function(inputdf){

     tissmaxTiss<-c()
     tissExpress<-c()

     for(i in c(1:nrow(inputdf))){
     
      tempmaxindex<-which(inputdf[i,] == max(inputdf[i,]))
      
      if(length(tempmaxindex >1)){
          tempmaxindex<- tempmaxindex[1]        
 } 
     tempmaxTiss<-colnames(inputdf)[tempmaxindex]
     tempmaxExpress<- inputdf[i,tempmaxindex]
     
      tissmaxTiss<-c(tissmaxTiss,tempmaxTiss)
      tissExpress<-c(tissExpress,tempmaxExpress)
           

}

return(list("Tissue"=tissmaxTiss,"Expression"=tissExpress))

}

#inputdf:avg expression matrix, rows are genes and cols are tissue


GeneGini<-function(input_df){
  genegini_all<-c()  
 for(i in c(1:nrow(input_df))){
         genegini<-Gini(input_df[i,])

   genegini_all<-c(genegini_all,genegini)


} 
 
return(genegini_all)

}
