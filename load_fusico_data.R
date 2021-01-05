##############################################################################/
##############################################################################/
#Data for the analysis and figure production of the fusicoccum project
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(rgdal)
library(plotrix)
library(classInt)
library(mapplots)
library(rgeos)
library(drc)
library(gdata)
library(RColorBrewer)

##############################################################################/
#loading the geographical data####
##############################################################################/

#the geographical layer used here were downloaded on the IGN (the French 
#Institut of Geographic and forest information) website: 
#http://professionnels.ign.fr/adminexpress and then turned into a .RData file
#The version of the dataset used here is the "Edition Novembre 2017"
#These data are under an open licence: 
#https://www.etalab.gouv.fr/wp-content/uploads/2014/05/Licence_Ouverte.pdf

#loading the different administrative unit levels in France
load("data/departe.RData")
load("data/regions.Rdata")


##############################################################################/
#loading the bioassay results data####
##############################################################################/

datafusamy<-read.table("data/fusicodat.txt",header=TRUE,sep="\t",
                       stringsAsFactors=TRUE)
datafuspop<-read.table("data/fusipopdat.txt",header=TRUE,sep="\t",
                       stringsAsFactors=TRUE)
datatemp<-read.table("data/fusico2020.txt",header=TRUE,sep="\t",
                     stringsAsFactors=TRUE)

##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/