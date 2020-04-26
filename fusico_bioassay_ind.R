##############################################################################/
##############################################################################/
#R code for analyzing bioassays of Fusicoccum amygdali isolates
##############################################################################/
##############################################################################/

#loading the data
source("load_fusico_data.R")


##############################################################################/
#Regression analysis of mycelial growth experiment for isolates####
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

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with the dose 
#range used with an absurdly high value
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-100
CompRez$STERR<-as.numeric(as.character(CompRez$STERR))
CompRez[is.na(CompRez$STERR),"STERR"]<-0


carbenfi<-CompRez[order(as.numeric(as.character(
  CompRez[CompRez$Subs_Act=="carbendazim",]$ED50))),]

plot(carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
     log="y",las=1,ylim=c(0.005,100),xlab="isolates",
     ylab="ED50 (mg/l)",col=carbenfi$popID,
     pch=c(19,19,19,22)[carbenfi$popID])

posi<-carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50+
  carbenfi[carbenfi$Subs_Act=="carbendazim",]$STERR
negi<-carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50-
  carbenfi[carbenfi$Subs_Act=="carbendazim",]$STERR

plotCI(c(1:63),
       carbenfi[carbenfi$Subs_Act=="carbendazim",]$ED50,
       ui=posi,
       li=negi,
       #ui=results$IC_up[results$sex=="female"],
       #li=results$IC_low[results$sex=="female"],
       add=TRUE,cex=0.1,pch=21,col=rgb(0,0,0,1),pt.bg=rgb(0.7,0.7,0.7,1),
       gap=0.00)


plot(sort(as.numeric(as.character(
  CompRez[CompRez$Subs_Act=="diethofencarb",]$ED50))),
  log="y",las=1,ylim=c(0.01,100))


##############################################################################/
#END
##############################################################################/