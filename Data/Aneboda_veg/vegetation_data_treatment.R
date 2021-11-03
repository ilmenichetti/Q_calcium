library(readxl)

Height<-read_excel("Aneboda_veg.xlsx", sheet=1)

Diameter<-read_excel("Aneboda_veg.xlsx", sheet=2)

Dead<-read_excel("Aneboda_veg.xlsx", sheet=3)


Damaged<-Dead!="*"
Damaged[is.na(Damaged)]=FALSE
Damaged<-as.data.frame(Damaged)
Damaged[,1:7]<-Dead[,1:7]

Damaged_years<-as.data.frame(cbind(Damaged$År, rowSums(Damaged[,8:23])))

aggregate(Damaged_years, by=list(Damaged_years$V1), FUN=sum)
