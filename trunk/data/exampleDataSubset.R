clinicalData<-read.table(file="exampleData.txt", sep="\t", quote="")
#Subset for asthma patients
exampleDataSubset<-clinicalData[which(clinicalData$Disease == "Asthma"),]
rm(clinicalData)
