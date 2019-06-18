###############################################################################
###############################################################################
#Code for the sampling map of Fusicoccum amygdali
###############################################################################
###############################################################################

#this code produce a map that is not displayed in the manuscript
source("load_fusico_data.R")


###############################################################################
#defining additionnal function for the mapping
###############################################################################


#function for a scale, found in "Auxiliary Cartographic Functions in R: 
#North Arrow, Scale Bar, and Label with a Leader Arrow", Tanimura et al 2007, 
#J of Statistical software
#The code has been slightly modified in order to convert the meter in km
scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
  if(missing(loc)) stop("loc is missing")
  if(missing(length)) stop("length is missing")
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- (x[c(1,3)]-loc[1])/1000
  labels <- append(labels,paste((x[5]-loc[1])/1000,unit))
  text(x[c(1,3,5)],y[4],labels=labels,adj=c(0.5,0),cex=division.cex)
}

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
  # checking arguments
  if(missing(loc)) stop("loc is missing")
  if(missing(size)) stop("size is missing")
  # default colors are white and black
  if(missing(cols)) cols <- rep(c("white","black"),8)
  # calculating coordinates of polygons
  radii <- rep(size/c(1,4,2,4),4)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
  # drawing polygons
  for (i in 1:15) {
    x1 <- c(x[i],x[i+1],loc[1])
    y1 <- c(y[i],y[i+1],loc[2])
    polygon(x1,y1,col=cols[i])
  }
  # drawing the last polygon
  polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                      (size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                      cex=cex)
}


###############################################################################
#Sampling site and resistance status of the sampled populations
###############################################################################

op<-par(mar=c(0,0,1,0))

#map summarizing the resistant and not resistant populations by department
temp<-datafuspop
colovec<-c(brewer.pal(9,"Blues")[6],brewer.pal(9,"Reds")[6])
#first we list the indices of the sampled department
ind_list<-which(departe@data$INSEE_DEP %in% 
                  colnames(table(temp$carbend_R,temp$departement)))
#because of strange departement denomination, we reorder the object
ind_list<-ind_list[c(3,4,1,2,5)]
#building the table of barycentre coordinates of the list of department
coorddep<-data.frame("longitude"=departe@polygons[ind_list[1]][[1]]@labpt[1],
                     "latitude"=departe@polygons[ind_list[1]][[1]]@labpt[2])
for (i in 2:length(ind_list)){
  coorddep<-rbind(coorddep, 
                  c("longitude"=departe@polygons[ind_list[i]][[1]]@labpt[1],
                    "latitude"=departe@polygons[ind_list[i]][[1]]@labpt[2]))
}
coorddep<-cbind(coorddep,"Res"=table(temp$carbend_R,temp$departement)[1,],
                "nonR"=if(dim(table(temp$carbend_R,temp$departement))[1]==1)
                  rep(0,dim(table(temp$carbend_R,temp$departement))[2])
                else table(temp$carbend_R,temp$departement)[2,],
                "nb_fields"=colSums(table(temp$carbend_R,temp$departement)))
plot(departe,border="grey60",lwd=0.1,main="")
plot(regions,add=TRUE,lwd=1.5)
draw.pie(x=coorddep$longitude,y=coorddep$latitude,
         z=cbind(coorddep$nonR,coorddep$Res),
         col=colovec,
         radius=(sqrt(coorddep$nb_fields)*18000),labels=NA,
         border="NA")
text(x=coorddep$longitude,y=coorddep$latitude,col="black",font=2,
     labels=as.character(coorddep$nb_fields),cex=1.5)

scalebar(c(191257.6,6080001),300000,"km",division.cex=1)

par(op)

#export pdf 6 x 6 inches


###############################################################################
#END
###############################################################################