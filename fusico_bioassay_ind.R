##############################################################################/
##############################################################################/
#R code for analyzing bioassays of Fusicoccum amygdali isolates
##############################################################################/
##############################################################################/

#loading the data
source("load_fusico_data.R")


##############################################################################/
#Regression analysis: carbendazim mycelial growth experiment for isolates####
##############################################################################/

datafusamyIND<-datafusamy[datafusamy$strain_type!="population",]
#first we extract the list of the different SA listed in the file
SAlist<-levels(datafusamyIND$active_substance)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=character(),STERR=character())
#we make a subselection of the data according to the SA
for (j in 1:length(SAlist)) {
  data_subSA<-datafusamyIND[datafusamyIND$active_substance==SAlist[j],]
  
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$perc_croiss>50,
                                  "strain_ID"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           ED50=character(),STERR=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "STERR"="unknown"))
  #we limit the dataset to the sample that reach somehow a IC of 50%
  SA.dat<-data_subSA[!(data_subSA$strain_ID %in% SA_rez),]
  SA.dat<-drop.levels(SA.dat)
  for (i in 1:dim(table(SA.dat$strain_ID))[1]) {
    temp.m1<-drm(perc_croiss~dose,
                 data=SA.dat[SA.dat$strain_ID==
                               names(table(SA.dat$strain_ID))[i],],
                 fct=LN.3())
    plot(temp.m1,ylim=c(0,120),xlim=c(0,150),
         main=paste(SAlist[j],names(table(SA.dat$strain_ID))[i]),
         col.main=j)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("Subs_Act"=SAlist[j],
                      "sample_ID"=names(table(SA.dat$strain_ID))[i],
                      "ED50"=as.character(temp[1]),
                      "STERR"=as.character(temp[2]))
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
  
}

#adding a column for the population the individuals were isolated from
CompRez<-data.frame(CompRez,"popID"=as.factor(substr(CompRez$sample_ID,1,7)))
#exporting the result as a text file
write.table(CompRez, file="output/results_fusicoIND.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Figure 3A: Distribution of the carbendazim ED50 on isolates####
##############################################################################/

#preparing the dataset, first we replace impossible value to compute by 
#the highest dose used in the bioassay
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-100
CompRez$STERR<-as.numeric(as.character(CompRez$STERR))
CompRez[is.na(CompRez$STERR),"STERR"]<-0

carbenfi<-CompRez[CompRez$Subs_Act=="carbendazim",]
carbenfi<-carbenfi[order(as.numeric(as.character(carbenfi$ED50))),]

#computing the values for the whiskers
posi<-carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50+
  carbenfi[carbenfi$Subs_Act=="carbendazim",]$STERR
negi<-carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50-
  carbenfi[carbenfi$Subs_Act=="carbendazim",]$STERR
#because we can't plot CI that reach negative values in a plot 
#with a log y-axes, we replace negative value by a very small value
negi[negi<0]<-0.001

#actual plotting
op<-par(mar=c(6.5,7,2,1),mfrow=c(2,1))
colov<-c("white","indianred1","black","dodgerblue")
plot(carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
     log="y",las=1,ylim=c(0.005,100),bty="n",axes=FALSE,
     ann=FALSE,col=colov[carbenfi$popID],
     bg=c("black","white","white","dodgerblue")[carbenfi$popID],
     pch=c(22,22,24,21)[carbenfi$popID],cex=1.5)
legend(45,35,col=colov,cex=1.5,x.intersp=0.5,y.intersp=0.4,
       pt.bg=c("black","white","white","dodgerblue"),
       legend=levels(carbenfi$popID),
       pch=c(22,22,24,21),bty="n")
box(lwd=2.5,lty=1)
axis(1,at=c(0,10,20,30,40,50,60),
     labels=c("0","10","20","30","40","50","60"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0.01,0.1,1,10,100),
     labels=c("0.01","0.1","1","10",">100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Isolates arranged by carbendazim\nED50 ascending order",
      ylab="ED50 (mg/l)",
      cex.lab=1.8,font.lab=2,line=5)
plotCI(c(1:63),
       carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
       ui=posi,
       li=negi,
       add=TRUE,cex=0.1,pch=21,col=rgb(0,0,0,1),pt.bg=rgb(0.7,0.7,0.7,1),
       gap=0.00)
points(carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
       pch=c(22,22,24,21)[carbenfi$popID],col=colov[carbenfi$popID],
       bg=c("black","white","white","dodgerblue")[carbenfi$popID],
       cex=1.5)
abline(h=50,lty=2,lwd=2)
text(-13,195,labels="(a)",cex=3,xpd=TRUE)

#export to .pdf 8 x 6 inches


##############################################################################/
#Figure 3B: Distribution of the diethofencarb ED50 on isolates####
##############################################################################/

#preparing the dataset, first we replace impossible value to compute by 
#the highest dose used in the bioassay
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-100
CompRez$STERR<-as.numeric(as.character(CompRez$STERR))
CompRez[is.na(CompRez$STERR),"STERR"]<-0

diethofe<-CompRez[CompRez$Subs_Act=="diethofencarb",]
diethofe<-diethofe[order(as.numeric(as.character(diethofe$ED50))),]

#computing the values for the whiskers
posi<-diethofe[diethofe$Subs_Act=="diethofencarb",]$ED50+
  diethofe[diethofe$Subs_Act=="diethofencarb",]$STERR
negi<-diethofe[diethofe$Subs_Act=="diethofencarb",]$ED50-
  diethofe[diethofe$Subs_Act=="diethofencarb",]$STERR
#because we can't plot CI that reach negative values in a plot 
#with a log y-axes, we replace negative value by a very small value
negi[negi<0]<-0.001

#actual plotting
colov<-c("white","indianred1","black","dodgerblue")
plot(diethofe[diethofe$Subs_Act=="diethofencarb",]$ED50,
     log="y",las=1,ylim=c(0.005,100),bty="n",axes=FALSE,
     ann=FALSE,col=colov[diethofe$popID],
     bg=c("black","white","white","dodgerblue")[diethofe$popID],
     pch=c(22,22,24,21)[diethofe$popID],cex=1.5)
box(lwd=2.5,lty=1)
axis(1,at=c(0,10,20,30,40,50,60),
     labels=c("0","10","20","30","40","50","60"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0.01,0.1,1,10,100),
     labels=c("0.01","0.1","1","10","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Isolates arranged by diethofencarb\nED50 ascending order",
      ylab="ED50 (mg/l)",
      cex.lab=1.8,font.lab=2,line=5)
plotCI(c(1:63),
       diethofe[diethofe$Subs_Act=="diethofencarb",]$ED50,
       ui=posi,li=negi,cex=0.1,pch=21,col=rgb(0,0,0,1),
       pt.bg=rgb(0.7,0.7,0.7,1),gap=0.00,add=TRUE)
points(diethofe[diethofe$Subs_Act=="diethofencarb",]$ED50,
       pch=c(22,22,24,21)[diethofe$popID],col=colov[diethofe$popID],
       bg=c("black","white","white","dodgerblue")[diethofe$popID],
       cex=1.5)
text(-13,195,labels="(b)",cex=3,xpd=TRUE)
par(op)

#export to .pdf 8 x 12 inches


##############################################################################/
#Figure SX: Correlation between ED50 carbendazim vs diethofencarb####
##############################################################################/

CarbVsDietho<-merge(carbenfi,diethofe,by="sample_ID")

op<-par(mar=c(6,6,2,1))
colov<-c("white","indianred1","black","dodgerblue")
plot(CarbVsDietho$ED50.y~CarbVsDietho$ED50.x,log="xy",bty="n",axes=FALSE,
     xlim=c(0.01,100),ylim=c(10,100),las=1,ann=FALSE,
     col=colov[CarbVsDietho$popID.x],
     bg=c("black","white","white","dodgerblue")[CarbVsDietho$popID.x],
     pch=c(22,22,24,21)[CarbVsDietho$popID.x],cex=1.5)
abline(v=50,lty=2,lwd=2)
legend(2,19,col=colov,cex=1.5,x.intersp=0.5,
       pt.bg=c("black","white","white","dodgerblue"),
       legend=levels(carbenfi$popID),
       pch=c(22,22,24,21),bty="o",box.col="transparent",bg="transparent")
box(lwd=2.5,lty=1)
axis(1,at=c(0.01,0.1,1,10,100),
     labels=c("0.01","0.1","1","10",">100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
axis(2,at=c(10,20,40,70,100),
     labels=c("10","20","40","70","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Carbendazim ED50 (mg/l)",
      ylab="Diethofencarb ED50 (mg/l)",
      cex.lab=1.8,font.lab=2,line=3.5)
par(op)
#text(0.0017,116,labels="C",cex=4,xpd=TRUE)

#export to .pdf 8 x 6 inches


##############################################################################/
#END
##############################################################################/



##############################################################################/
#Regression analysis: carbendazim mycelial growth experiment for isolates####
##############################################################################/

datafusamyIND<-datatemp[datatemp$strain_type!="population",]
#first we extract the list of the different SA listed in the file
SAlist<-levels(datafusamyIND$active_substance)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=character(),STERR=character())
#we make a subselection of the data according to the SA
for (j in 1:length(SAlist)) {
  data_subSA<-datafusamyIND[datafusamyIND$active_substance==SAlist[j],]
  
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$perc_croiss>50,
                                  "strain_ID"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           ED50=character(),STERR=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "STERR"="unknown"))
  #we limit the dataset to the sample that reach somehow a IC of 50%
  SA.dat<-data_subSA[!(data_subSA$strain_ID %in% SA_rez),]
  SA.dat<-drop.levels(SA.dat)
  for (i in 1:dim(table(SA.dat$strain_ID))[1]) {
    temp.m1<-drm(perc_croiss~dose,
                 data=SA.dat[SA.dat$strain_ID==
                               names(table(SA.dat$strain_ID))[i],],
                 fct=LN.3())
    plot(temp.m1,ylim=c(0,120),xlim=c(0,450),
         main=paste(SAlist[j],names(table(SA.dat$strain_ID))[i]),
         col.main=j)
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("Subs_Act"=SAlist[j],
                      "sample_ID"=names(table(SA.dat$strain_ID))[i],
                      "ED50"=as.character(temp[1]),
                      "STERR"=as.character(temp[2]))
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
  
}

#adding a column for the population the individuals were isolated from
CompRez<-data.frame(CompRez,"popID"=as.factor(substr(CompRez$sample_ID,1,7)))
#exporting the result as a text file
write.table(CompRez, file="output/results_fusitemp.txt",
            sep="\t",quote=FALSE,row.names=FALSE)



