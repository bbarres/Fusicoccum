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
#Figure 2: Distribution of the carbendazim ED50 on isolates (myc. growth)####
##############################################################################/

#preparing the dataset, first we replace impossible value to compute by 
#the highest dose used in the bioassay
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-100
CompRez$STERR<-as.numeric(as.character(CompRez$STERR))
CompRez[is.na(CompRez$STERR),"STERR"]<-0

carbenfi<-CompRez[order(as.numeric(as.character(
  CompRez[CompRez$Subs_Act=="carbendazim",]$ED50))),]
#computing the values for the whiskers
posi<-carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50+
  carbenfi[carbenfi$Subs_Act=="carbendazim",]$STERR
negi<-carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50-
  carbenfi[carbenfi$Subs_Act=="carbendazim",]$STERR
#because we can't plot CI that reach negative values in a plot 
#with a log y-axes, we replace negative value by a very small value
negi[negi<0]<-0.001

#actual plotting
op<-par(mar=c(6,6,2,1))
colov<-c("black","indianred1","black","dodgerblue")
plot(carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
     log="y",las=1,ylim=c(0.005,100),bty="n",axes=FALSE,
     ann=FALSE,col=colov[carbenfi$popID],
     bg=c("black","white","white","dodgerblue")[carbenfi$popID],
     pch=c(22,22,24,21)[carbenfi$popID],cex=1.5)
legend(1,130,col=colov,cex=1.5,x.intersp=0.5,
       legend=levels(carbenfi$popID),
       pch=c(15,22,2,19),bty="n")
box(lwd=2.5,lty=1)
axis(1,at=c(0,10,20,30,40,50,60),
     labels=c("0","10","20","30","40","50","60"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0.01,0.1,1,10,100),
     labels=c("0.01","0.1","1","10","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Isolates arranged by ED50 ascending order",
      ylab="ED50 (mg/l)",
      cex.lab=1.8,font.lab=2,line=3.5)
plotCI(c(1:63),
       carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
       ui=posi,
       li=negi,
       add=TRUE,cex=0.1,pch=21,col=rgb(0,0,0,1),pt.bg=rgb(0.7,0.7,0.7,1),
       gap=0.00)
points(carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
       pch=c(22,22,24,21)[carbenfi$popID],col=colov[carbenfi$popID],
       bg=c("black","white","white","dodgerblue")[carbenfi$popID],
       add=TRUE,cex=1.5)
par(op)




plot(sort(as.numeric(as.character(
  CompRez[CompRez$Subs_Act=="diethofencarb",]$ED50))),
  log="y",las=1,ylim=c(0.01,100))


##############################################################################/
#END
##############################################################################/