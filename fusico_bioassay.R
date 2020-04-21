##############################################################################/
##############################################################################/
#R code for analyzing the ouput of bioassays of Fusicoccum amygdali
##############################################################################/
##############################################################################/

#loading the data
source("load_fusico_data.R")


##############################################################################/
#Regression analysis of percentage of sporulation for the populations####
##############################################################################/

datafusamyPOP<-datafusamy[datafusamy$strain_type=="population",]
#first we extract the list of the different SA listed in the file
SAlist<-levels(datafusamyPOP$active_substance)
CompRezPOP<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=character(),STERR=character())
recap.mod<-list()
#we make a subselection of the data according to the SA
for (j in 1:length(SAlist)) {
  data_subSA<-datafusamyPOP[datafusamyPOP$active_substance==SAlist[j],]
  
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$perc_croiss>50,
                                  "strain_ID"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           ED50=character(),STERR=character(),
                           ED99.9=character(),STERR=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "STERR"="unknown",
                           "ED99.9"=paste(">",max(data_subSA$dose),sep=""),
                           "STERR"="unknown"))
  #we limit the dataset to the sample that reach somehow a IC of 50%
  SA.dat<-data_subSA[!(data_subSA$strain_ID %in% SA_rez),]
  SA.dat<-drop.levels(SA.dat)
  for (i in 1:dim(table(SA.dat$strain_ID))[1]) {
    temp.m1<-drm(perc_croiss~dose,
                 data=SA.dat[SA.dat$strain_ID==
                               names(table(SA.dat$strain_ID))[i],],
                 fct=LN.3())
    plot(temp.m1,main=names(table(SA.dat$strain_ID))[i],
         ylim=c(0,120),xlim=c(0,150))
    temp<-ED(temp.m1,c(50,0.1),type="absolute")
    tempx<-data.frame("Subs_Act"=SAlist[j],
                      "sample_ID"=names(table(SA.dat$strain_ID))[i],
                      "ED50"=as.character(temp[1,1]),
                      "STERR"=as.character(temp[1,2]),
                      "ED99.9"=as.character(temp[2,1]),
                      "STERR"=as.character(temp[2,2]))
    REZSA<-rbind(REZSA,tempx)
    recap.mod<-c(recap.mod,temp.m1)
  }
  CompRezPOP<-rbind(CompRezPOP,REZSA)
  
}

#exporting the result as a text file
write.table(CompRezPOP, file="output/results_fusicoPOP.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


datafusamyPOP<-datafusamy[datafusamy$strain_type=="population",]
data_subSA<-datafusamyPOP[datafusamyPOP$active_substance=="carbendazim",]
#simplier way of doing approximately the same thing
temp.m1<-drm(perc_croiss~dose,
             data=data_subSA,curveid=strain_ID,
             fct=LN.3())
plot(temp.m1,xlim=c(0,150),lwd=2)
summary(temp.m1)
modelFit(temp.m1)
compParm(temp.m1,"e")
temp<-ED(temp.m1,c(50,0.1),type="absolute")


#just a small graphic to gain insight on the first round of results first, 
#we replace the ED50 that were too high to be evaluated with the dose range 
#used with an absurdly high value
CompRezPOP$ED50<-as.numeric(as.character(CompRezPOP$ED50))
CompRezPOP[is.na(CompRezPOP$ED50),"ED50"]<-100
barplot(as.numeric(as.character(CompRezPOP$ED50)),
        ylim=c(0,50),col=CompRezPOP$Subs_Act)

plot(sort(as.numeric(as.character(
  CompRezPOP[CompRezPOP$Subs_Act=="carbendazim",]$ED50))),
  ylog=TRUE)


##############################################################################/
#Regression analysis of mycelial growth experiment####
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

#exporting the result as a text file
write.table(CompRez, file="output/results_fusicoIND.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with the dose 
#range used with an absurdly high value
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-100
plot(sort(as.numeric(as.character(
  CompRez[CompRez$Subs_Act=="carbendazim",]$ED50))),
  log="y",las=1,ylim=c(0.01,100))
plot(sort(as.numeric(as.character(
  CompRez[CompRez$Subs_Act=="diethofencarb",]$ED50))),
  log="y",las=1,ylim=c(0.01,100))


##############################################################################/
#END
##############################################################################/