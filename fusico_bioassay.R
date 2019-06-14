###############################################################################
###############################################################################
#R code for analyzing the ouput of bioassays of Fusicoccum amygdali
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(drc)
library(plotrix)
library(gdata)

#loading the data
datafusamy<-read.table("data/fusicodat.txt",header=TRUE,sep="\t")


###############################################################################
#Regression analysis of mycelial growth experiment
###############################################################################

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=character())
#we make a subselection of the data according to the SA
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  SA_rez<-as.character(data_subSA[data_subSA$dose==max(data_subSA$dose) 
                                  & data_subSA$rslt_03>50,
                                  "ech_id"])
  ifelse(length(SA_rez)==0,
         REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                           ED50=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50"=paste(">",max(data_subSA$dose),sep="")))
  #we limit the dataset to the sample that reach somehow a IC of 50%
  SA.dat<-data_subSA[!(data_subSA$ech_id %in% SA_rez),]
  SA.dat<-drop.levels(SA.dat)
  for (i in 1:dim(table(SA.dat$ech_id))[1]) {
    temp.m1<-drm(rslt_03~dose,
                 data=SA.dat[SA.dat$ech_id==names(table(SA.dat$ech_id))[i],],
                 fct=LL.4())
    plot(temp.m1,main=names(table(SA.dat$ech_id))[i],
         ylim=c(0,110),xlim=c(0,50))
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("Subs_Act"=SAlist[j],
                      "sample_ID"=names(table(SA.dat$ech_id))[i],
                      "ED50"=as.character(temp[1]))
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
  
}

#exporting the result as a text file
write.table(CompRez, file="output/results_cerco.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with the dose 
#range used with an absurdly high value
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))
CompRez[is.na(CompRez$ED50),"ED50"]<-10000
barplot(as.numeric(as.character(CompRez$ED50)),
        ylim=c(0,50),col=CompRez$Subs_Act)


###############################################################################
#END
###############################################################################