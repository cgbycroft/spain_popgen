#################################################
# Download map data and plot points on maps
# coloured by some arbitrary quantitative value.
#
# Author: Clare Bycroft
# Date: February 2015
#################################################
# Version requirements:
# 
# R: 3.1.0 or higher
# Packages (see below)
#################################################

#################################################
#  Required packages.
#################################################

packages=c("colortools_0.1.5","rgdal_0.8-16","colorspace_1.2-4","latticeExtra_0.6-26",
           "RColorBrewer_1.0-5","lattice_0.20-29","maptools_0.8-30","maps_2.3-7",
           "raster_2.2-31","sp_1.0-15","rgeos_0.3-5")

#  If you don't have  with :
#  install.packages(pkgs=packages)

library(rgeos) ### rgeos_0.3-5
library(raster) ### raster_2.2-31
library(maps) ### maps_2.3-7
library(maptools) ### maptools_0.8-30
library(sp) ### sp_1.0-15
library(lattice) ### lattice_0.20-29
library(latticeExtra) ### latticeExtra_0.6-26
library(colorspace) ### colorspace_1.2-4
library(grid) ### (R base package: 3.1.0 or higher)
library(rgdal) ### rgdal_0.8-16
library(colortools) ### colortools_0.1.5

#################################################
# some helpful functions
#################################################

colourScale <- function(x,colourSet= c("green","yellow","red"),nBreaks=200,fixedLims=NULL){
                Colors <- colorRampPalette(colourSet)(nBreaks)
                if (!is.null(fixedLims)) scale <- seq(fixedLims[1],fixedLims[2]+0.0001,length.out=nBreaks) else scale <- seq(min(x),max(x)+0.0001,length.out=nBreaks)
                cols <- c()
                for (val in x) cols <- c(cols,Colors[which(scale>val)[1]])
return(cols)
}

makeScaleWithHist <- function(values,colourSet,scalenum=100,nbreaks=70,fixedLims=NULL,binWidth=NULL,xlims=NULL,bgColour="transparent"){
    par(mar=c(0,0,0,0))
    if(is.null(fixedLims)) colindex<-matrix(seq(min(values),max(values),length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
    if(!is.null(fixedLims)) colindex<-matrix(seq(fixedLims[1],fixedLims[2],length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
   
    if(!is.null(binWidth)) getBreaks <- seq(min(values),max(values)+binWidth,by=binWidth) else getBreaks <- seq(min(values),max(values),length.out=nbreaks)
   
    h <- hist(values,breaks=getBreaks,plot=F)
    histMax <- max(h$density) + max(h$density)/20
    if(is.null(xlims)) xlims=c(min(values),max(values))
    
    hist(values,breaks=getBreaks,xlim=xlims,axes=F,xlab=NULL,ylab=NULL,main=NULL,probability=T,border="gray")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bgColour,border=NA)
    #axis(1,cex.axis=0.5,tck=-0.1)
    
    colours <- colourScale(colindex[,1],colourSet)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=length(colindex))
    #scalephysicalpos <- seq(1,nbreaks,length.out=length(colindex))
    scalephysicalpos <- seq(1,nbreaks,length.out=length(colindex))
   # image(1:nbreaks,1,colindex,xaxt="n",xlab="",yaxt="n",ylab="",col=colours,zlim=range(colindex))
     image(scalelocs,c(-max(h$density)/20,histMax),colindex,xaxt="n",xlab="",yaxt="n",ylab="",col=colours,zlim=range(colindex),axes=F,add=T)     
    # axis(1,las=1,cex.axis=0.6,tcl=-0.2,padj=-3)
    # axis(1,las=1,cex.axis=0.6)
    
    
    par(mgp=c(0,0,0))
    if(is.null(fixedLims)) hist(values,breaks=getBreaks,col=My.add.alpha("white",0.7),axes=F,xlab=NULL,ylab=NULL,main=NULL,border="transparent",probability=T,add=T,xlim=xlims)
    if(!is.null(fixedLims)) hist(values,breaks=getBreaks,col=My.add.alpha("white",0.8),
                                 axes=F,xlab=NULL,ylab=NULL,main=NULL,border="transparent",xlim=c(fixedLims[1],fixedLims[2]),probability=T,add=T)
    axis(1,cex.axis=0.5,tck=-0.1,at=signif(seq(xlims[1],xlims[2],length.out=7),2))
     #box(lwd=0.5)    
    par(mar=opar$mar,mgp=opar$mgp)
}


#################################################
# download country maps from here: http://www.gadm.org/country
#################################################


setwd("~/Documents/ClareDPhil/DPhil/Spain/mapping")

load('~/Downloads/CHN_adm3.Rdata',verbose=T)

chinaMap <- gadm

# Combine some polygons using gUnaryunion

chinaMapName1 <- gUnaryUnion(chinaMap,chinaMap$NAME_1)

# This creates a SpatialPolygons object, so to map it you need to create a SpatialPolygonsDataFrame
npolygons <- length(chinaMapName1@polygons)
Data <- as.data.frame(rep(1,npolygons))
  rownames(Data) <- sapply(chinaMapName1@polygons,FUN=slot,name="ID")
  colnames(Data) <- "dummy"

chinaMapName1df <- SpatialPolygonsDataFrame(chinaMapName1,data=Data,match.ID=TRUE)

# Plot the map simply

spplot(chinaMapName1df,zcol="dummy",regioncol="gray",colorkey=F)

# Colour polygons by some measure <measure> which is recorded as a column in chinaMapName1df@data
chinaMapName1df@data$measure <- runif(npolygons)
spplot(chinaMapName1df[1:5,],zcol="measure",regioncol=c("red","blue"),colorkey=T)

                      
# It takes a while to plot, so do we need all of that resolution? (7,000 points describe one polygon!)
      # Use gSimplify
                      
chinaMapName1simple <- gSimplify(chinaMapName1df,tol=0.1,topologyPreserve=T)
chinaMapName1simpledf  <- SpatialPolygonsDataFrame(chinaMapName1simple,data=chinaMapName1df@data,match.ID=T)
colourPal <- colorRampPalette(colors=c("red","blue"))(50)     
                      # look at colours using pizza!!!
                      pizza(colourPal)
                      
chinaMapRegions <- spplot(chinaMapName1simpledf,zcol="measure",col.regions=colourPal,colorkey=T,col="gray")

# What about the rest of the world?

     # What are the countries called ?
country_names <- getData('ISO3')
                      
### Download country data (only do from ****** to ****** once)

#******                     
system('mkdir ./world_map_files')                      
getData('countries',path='./world_map_files')
countries$dummy <- 0
proj4string(countries) <- CRS(proj4string(countries))

for (column in colnames(countries@data)){
  if (!class(countries@data[,column])%in%c("character","factor")) next
  print(column)
  Encoding(levels(countries@data[,column])) <- "latin1"
}
                      
#******
                      
load('./world_map_files/countries.Rdata')                      

# Just get polygons for Asia
Asia <- countries[countries@data$CONTINENT=="Asia",]                      
spplot(Asia,zcol="dummy")

# This is a bit too big, so get just countries around China
                      
subbox <- countries[(countries@data$COUNTRY=="China")&(countries@data$CONTINENT=="Asia"),]@bbox

subboxPolygon <- SpatialPolygons(list(Polygons(list(Polygon(rbind(subbox[,1],c(subbox[1,1],subbox[2,2]),subbox[,2],c(subbox[1,2],subbox[2,1]),subbox[,1]))),"subbox")),
                                 proj4string=CRS(proj4string(countries)))

oversub <- which(!is.na(over(countries,subboxPolygon)))
countriesOverChina <- countries[oversub,]

countriesOverChinasimple <- gSimplify(countriesOverChina,tol=0.1,topologyPreserve=T)
countriesOverChinasimpledf  <- SpatialPolygonsDataFrame(countriesOverChinasimple,data=countriesOverChina@data,match.ID=T)

chinaMapCountry <- spplot(countriesOverChinasimpledf,zcol="dummy",xlim=c(subbox[1,1],subbox[1,2]),ylim=c(subbox[2,1],subbox[2,2]),
                        col.regions="transparent",col="black",
                        lwd=1.2,colorkey=F)
          # col = colour of borders
          # col.regions = colour inside polygons


# Put map of China regions on top of context map
colouredPlot <- chinaMapRegions + chinaMapCountry

png('ChinaColouredRandom.png',height=2000,width=2000,res=150)
  print(colouredPlot)
dev.off()

# A basic plot for putting points on (no colour key)

BasicPlot <- chinaMapCountry + update(chinaMapRegions,col.regions="lightgray")

png('ChinaBasicPlot.png',height=2000,width=2000,res=150)
  print(BasicPlot)
dev.off()

# How to add points

# Generate some random points
X <- runif(100,subbox[1,1],subbox[1,2])
Y <- runif(100,subbox[2,1],subbox[2,2])
coords <- cbind(X,Y)

# Create SpatialPoints object

randPoints <- SpatialPoints(coords,proj4string=CRS(proj4string(chinaMapName1simpledf)))
someNumVector <- coords[,1]

cols <- colourScale(someNumVector,colourSet=c("red","blue"),nBreaks=50)

mappoints <- list("sp.points",randPoints,cex=2,pch=16,col=cols)
#maplines <- list("sp.lines",...)

myPlot <- chinaMapCountry + update(chinaMapRegions,col.regions="lightgray",sp.layout=list(mappoints))


# Plot scale with map
mapViewport = viewport(0.5,1,just=c("centre","top"))
scaleViewport = viewport(0.05,0.03,height=0.05,width=unit(0.9,"npc"),just=c("left","bottom"))
  
png('ChinapPointsPlot.png',height=2000,width=2000,res=150)
plot.new()
grid.newpage()
    pushViewport(mapViewport)
      print(myPlot,newpage=F)
    upViewport(1)
    pushViewport(scaleViewport)
      par(new=T,fig=gridFIG())
      makeScaleWithHist(values=someNumVector,colourSet=c("red","blue"),binWidth=1,scalenum=50)
      
dev.off()

#### get some coordinates and jitter

# get coordinates for towns (one for each polygon)
centres <- t(sapply(chinaMap@polygons,FUN=slot,name="labpt"))

# e.g Shanghai
shanghaiCoord <- centres[grep("Shanghai",chinaMap@data$NAME_3),]

# sample points within Shanhai
randShangPoints <- spsample(chinaMap[grep("Shanghai",chinaMap@data$NAME_3),], n=100, type='random')

# plot on map
mappoints <- list("sp.points",randShangPoints,cex=1,pch=16,col=cols)
myPlot <- chinaMapCountry + update(chinaMapRegions,col.regions="lightgray",sp.layout=list(mappoints))

png('ChinaShanghaiPointsPlot.png',height=2000,width=2000,res=150)
  print(myPlot)
dev.off()

# sample points within some province
randGanPoints <- spsample(chinaMap[grep("Gansu",chinaMap@data$NAME_1),], n=100, type='random')

# plot on map
mappoints <- list("sp.points",randGanPoints,cex=1,pch=16,col=cols)
myPlot <- chinaMapCountry + update(chinaMapRegions,col.regions="lightgray",sp.layout=list(mappoints))

png('ChinaGanPointsPlot.png',height=2000,width=2000,res=150)
  print(myPlot)
dev.off()
