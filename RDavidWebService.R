##Tutorial for using RDAVIDWebService 

##Helpful Site https://gist.github.com/svigneau/9699239

##You first have to register at https://david.ncifcrf.gov/webservice/register.htm




library(rJava)
library(RDAVIDWebService)

##Create david object 

david<- DAVIDWebService(email="yankee@uchc.edu",
	url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")






