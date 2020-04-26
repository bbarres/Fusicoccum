##############################################################################/
##############################################################################/
#R code for analyzing the bioassays of Fusicoccum amygdali populations
##############################################################################/
##############################################################################/

#loading the data
source("load_fusico_data.R")


##############################################################################/
#Regression analysis of percentage of sporulation for the populations####
##############################################################################/

#data selection and prepartion
datafusamyPOP<-datafusamy[datafusamy$strain_type=="population",]
datafusamyPOP<-drop.levels(datafusamyPOP)

#simplier way of doing approximately the same thing
temp.m1<-drm(perc_croiss~dose,
             data=datafusamyPOP,curveid=strain_ID,
             fct=LN.4())
summary(temp.m1)
compParm(temp.m1,"e")
result_pop<-ED(temp.m1,c(50),type="absolute")


##############################################################################/
#Figure 1: regression curves for the ####
##############################################################################/

plot(temp.m1,xlim=c(0,150),lwd=2,col=c(1,1,1,1,2,1,1,1))


#export to .pdf 6 x 6 inches


##############################################################################/
#alternatively, we can analyze the population one by one####
##############################################################################/

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
                           ED50=character(),STERR=character()),
         REZSA<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=SA_rez,
                           "ED50"=paste(">",max(data_subSA$dose),sep=""),
                           "STERR"="unknown")
  )
  #we limit the dataset to the sample that reach somehow a IC of 50%
  SA.dat<-data_subSA[!(data_subSA$strain_ID %in% SA_rez),]
  SA.dat<-drop.levels(SA.dat)
  for (i in 1:dim(table(SA.dat$strain_ID))[1]) {
    temp.m1<-drm(perc_croiss~dose,
                 data=SA.dat[SA.dat$strain_ID==
                               names(table(SA.dat$strain_ID))[i],],
                 fct=LN.4())
    plot(temp.m1,main=names(table(SA.dat$strain_ID))[i],
         ylim=c(0,120),xlim=c(0,150),type="confidence")
    temp<-ED(temp.m1,c(50),type="absolute")
    tempx<-data.frame("Subs_Act"=SAlist[j],
                      "sample_ID"=names(table(SA.dat$strain_ID))[i],
                      "ED50"=as.character(temp[1,1]),
                      "STERR"=as.character(temp[1,2]))
    REZSA<-rbind(REZSA,tempx)
    recap.mod[[paste(names(table(SA.dat$strain_ID))[i])]]<-temp.m1
  }
  CompRezPOP<-rbind(CompRezPOP,REZSA)
  
}

#exporting the result as a text file
write.table(CompRezPOP, file="output/results_fusicoPOP.txt",
            sep="\t",quote=FALSE,row.names=FALSE)

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
#END
##############################################################################/