library(ineq)
tissue_colval<-function(input){
  
  return(which(input == as.character(tm$smts)))
  
}


Tissue_Average<-function(input){
  tissavg<-c()
  
  for( i in c(1:nrow(pt5))){
     tempmean<-round(mean(pt5[i,input]),2)
   

    tissavg<-c(tissavg,tempmean)
   
}
  

 return(tissavg)
  
}

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

GeneGini<-function(input_df){
  genegini_all<-c()  
 for(i in c(1:nrow(input_df))){
         genegini<-Gini(input_df[i,])

   genegini_all<-c(genegini_all,genegini)


} 
 
return(genegini_all)

}
