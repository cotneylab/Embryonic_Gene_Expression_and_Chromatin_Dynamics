cout<-0
cin<-0
for (i in dss_filt){
  
  if (i >= 2000){
    cout<-cout+1
  }
  else if (i <2000){
    cin<-cin+1
  }
  
}
slices<-c(cout,cin)
lbls<-c("TSS distal","TSS proximal")
pie(slices,labels=lbls,main="dat1")