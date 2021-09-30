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

#data selection and preparation
datafusamyPOP<-datafusamy[datafusamy$strain_type=="population",]
datafusamyPOP<-drop.levels(datafusamyPOP)

#simplier way of doing approximately the same thing
temp.m1<-drm(perc_croiss~dose,
             data=datafusamyPOP,curveid=strain_ID,
             fct=LN.4())
summary(temp.m1)
compParm(temp.m1,"e")
result_pop<-ED(temp.m1,c(50),type="absolute")
result_pop<-data.frame("pop_ID"=sort(levels(datafusamyPOP$strain_ID)),
                       result_pop)

#exporting the result as a text file
write.table(result_pop, file="output/results_fusicoPOP.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#Figure 2: regression curves for the population bioassays with carbendazim####
##############################################################################/

op<-par(mar=c(6,7,2,1))
colov<-c("black","indianred1")
plot(temp.m1,xlim=c(0,200),lwd=1.7,col=colov[c(1,1,1,1,2,1,1,1)],
     pch=c(15:17,21:25),lty=c(1:8),
     bty="n",axes=FALSE,ann=FALSE,legend=FALSE)
# plot(temp.m1,xlim=c(0,200),lwd=1.7,col=colov[c(1,1,1,1,2,1,1,1)],
#      pch=c(15:17,21:25),lty=c(1:8),type="bars",
#      bty="n",axes=FALSE,ann=FALSE,legend=FALSE,add=TRUE)
legend(5,133,col=colov[c(1,1,1,1,2,1,1,1)],
       legend=levels(as.factor(temp.m1$parNames[[3]])),
       lty=c(1:8),lwd=1.7,pch=c(15:17,21:25),bty="n")
box(lwd=2.5,lty=1)
axis(1,at=c(0.001,0.01,0.1,1,10,100),
     labels=c("0","0.01","0.1","1","10","100"),
     cex.axis=1.1,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,20,40,60,80,100,120),
     labels=c("0","20","40","60","80","100","120"),
     cex.axis=1.1,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Dose (mg/L)",ylab="Relative % of germ tube length",
      cex.lab=1.5,font.lab=2,line=4)
par(op)

#export to .pdf 6.5 x 6 inches


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
write.table(CompRezPOP, file="output/results_fusicoPOP2.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#END
##############################################################################/