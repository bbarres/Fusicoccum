###############################################################################
###############################################################################
#Data for the analysis and figure production of the fusicoccum project
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(rgdal)
library(plotrix)
library(classInt)
library(mapplots)
library(rgeos)
library(drc)
library(gdata)


###############################################################################
#loading the geographical data
###############################################################################

#the geographical layer used here were downloaded on the IGN (the French 
#Institut of Geographic and forest information) website: 
#http://professionnels.ign.fr/adminexpress and then turned into a .RData file
#The version of the dataset used here is the "Edition Novembre 2017"
#These data are under an open licence: 
#https://www.etalab.gouv.fr/wp-content/uploads/2014/05/Licence_Ouverte.pdf

#loading the different administrative unit levels in France
load("data/departe.RData")
load("data/regions.Rdata")

#isolate the information in the spatial data on the communes
db_commu<-commu@data
summary(db_commu)


###############################################################################
#loading the bioassay results data
###############################################################################

datafusamy<-read.table("data/fusicodat.txt",header=TRUE,sep="\t")
datafuspop<-read.table("data/fusipopdat.txt",header=TRUE,sep="\t")


###############################################################################
#END
###############################################################################