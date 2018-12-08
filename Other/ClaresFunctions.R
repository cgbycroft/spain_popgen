##########################################
#  Various functions used in plotting etc.
##########################################

packages=c("gplots",
"ggplot2",
"scales",
"RColorBrewer",
"Matrix",
"stringr",
"rgeos",
"rgdal",
"maps",
"maptools",
"sp",
"lattice",
"latticeExtra",
"gridExtra",
"colorspace",
"grid",
"XML",
"ape",
"geiger",
"phangorn",
"dendroextras",
"colortools",
"plyr",
"ggdendro",
"plotrix",
"gridBase",
"nnet",
"abind",
"aberrant")

#install.packages(pkgs=packages)
library(gplots)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(Matrix)
library(stringr)
library(rgeos)
library(rgdal)
library(gdalUtils)
library(insol)
library(maps)
library(maptools)
library(sp)
library(colorspace)
library(grid)
library(XML)
library(ape)
library(phytools)
library(geiger)
library(phangorn)
library(dendroextras)
library(dendextend)
library(dendextendRcpp) # a much faster version of dendextend!
library(colortools)
library(plyr)
library(ggdendro)
library(plotrix)
library(gridBase)
library(nnet)
library(abind)
library(aberrant)
library(fBasics)
library(rasterVis)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(sfsmisc)
#library(dplyr)
library(RColorBrewer)
library(RgoogleMaps)
library(dplyr)
library(ade4)
library(igraph)
library(qgraph)
library(xtable)
library(viridis)

source("~/Documents/ClareDPhil/DPhil/Spain/finestructureLibrary/FinestructureLibrary.R")

# Note - only works for Spanish data as it uses legend.master
load('~/Documents/ClareDPhil/DPhil/Spain/samples/legend.master.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/samples/legend.master.hospitals.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/samples/legend.master.context.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/samples/SPAIN.A.attributes.all.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/samples/SPAIN.B.attributes.all.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/samples/SP_PR_HM_NA.attributes.all.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/samples/numToColourForGaussian_v1.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/espComMap.Rdata') 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/espMapTotal.Rdata') 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/espProvTotal.Rdata') 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/espMuniTotal.Rdata') 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/espProvTotal2.Rdata') 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/espMuniTotal2.Rdata') 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/ALL.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/ALL_context.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/ALL_context_map.Rdata')
load("~/Documents/ClareDPhil/DPhil/Spain/mapping/SPAIN.A2.lines.tight.Rdata") 
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/point.concordance.A.Rdata')
load('~/Documents/ClareDPhil/DPhil/Spain/mapping/point.concordance.PR_HM_NA.Rdata')

#####################
# Set some global variables
#####################
spainDIR="~/Documents/ClareDPhil/DPhil/Spain"
banksia <- 'banksia.well.ox.ac.uk:/well/donnelly/spain_structure'
eucalyptus <- 'eucalyptus.well.ox.ac.uk:/well/donnelly/spain_structure'
coolibah <- 'coolibah.well.ox.ac.uk:/well/donnelly/spain_structure'

opar <- par()

#####################
# Some generic useful functions
#####################

# append lists
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

#Extract Legend 
g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

# scp to and from banksia, or other servers

SCP <- function(filename,direction="from",recievedir="",spain=T,verbose=F,server=eucalyptus){
  if(direction=="from") {
    if (spain==T) CopyCommand <- paste(server,"/",filename," ",spainDIR,"/",recievedir,sep="")
    if (spain==F) CopyCommand <- paste(filename," ",recievedir,sep="")
  }
  if(direction=="to") {
    if (spain==T) CopyCommand <- paste(spainDIR,"/",filename," ",server,"/",recievedir,sep="")
    if (spain==F) CopyCommand <- paste(filename," ",recievedir,sep="")
  }
  if(verbose==T) CopyCommand <- paste('scp -v ',CopyCommand,sep="") else CopyCommand <- paste('scp ',CopyCommand,sep="") 
  print(CopyCommand)
  system(CopyCommand)
}

# add grayness to colours

My.add.alpha <- function(col, alpha=1, valueChange=1,saturationChange=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb), 2, 
                     function(x) {
                       Hsv=rgb2hsv(x)
                       hsv(Hsv[1],Hsv[2]*saturationChange,Hsv[3]*valueChange,alpha=alpha)
                     })
}

My.add.alpha2 <- function(col, alpha=1, value=1,saturation=1){
  if(missing(col))
    # saturation and value are absolute values for the output (rather than fold changes)
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb), 2, 
                     function(x) {
                       Hsv=rgb2hsv(x)
                       hsv(Hsv[1],saturation,value,alpha=alpha)
                     })
}

# get colour scales from numeric vector
colourScale <- function(x,colourSet= c("green","yellow","red"),nBreaks=200,fixedLims=NULL,colSpace="Lab",...){
                if(class(colourSet)=="function") Colors = colourSet(nBreaks) else Colors <- colorRampPalette(colourSet,space=colSpace,...)(nBreaks)
                if (!is.null(fixedLims)) scale <- seq(fixedLims[1],fixedLims[2]+0.0001,length.out=nBreaks) else scale <- seq(min(x),max(x)+0.0001,length.out=nBreaks)
                cols <- c()
                for (val in x) cols <- c(cols,Colors[which(scale>val)[1]])
return(cols)
}

# List of colourblind-friendly colours from http://mkweb.bcgsc.ca/colorblind/img/colorblindness.palettes.trivial.png
getColBlindCol <- function(want=NULL){
    o=rgb2hsv(230,159,0) # orange
    s=rgb2hsv(86,180,233) # sky blue
    g=rgb2hsv(0,158,115) # bluegreen
    y=rgb2hsv(240,228,66) # yellow
    b=rgb2hsv(0,114,178) # blue
    v=rgb2hsv(213,94,0) # vermillion
    r=rgb2hsv(204,121,167) # red/pink
    if(is.null(want)) want = c("y","o","v","r","b","s","g")
    cols = sapply(want,function(x) {
        m=get(x)
        return(hsv(m[1],m[2],m[3]))
        })
    return(cols)
}

primaryColors = c("red","blue","green","darkgreen","yellow","cyan","orange","purple","maroon4","darkgreen","darkblue","skyblue","orange3")



# rotate matrix 90 degrees

rotate90=function(m){
  mt=t(m)
  out=mt[,ncol(mt):1]  
}

# aggregate rows of a matrix over group categories using an arbitrary function

agg <- function(data,byGroup,myFunction="sum",useNames=T,...){
    # myFunction is a name of a function, or an actual function, with the first argument being a vector (i.e a subset of data). Default function is sum.
    # ... are further arguments to be fed into myFunction via do.call()
    # byGroup is a vector of categories, in the order of the columns of data
    uniqueGroups = unique(byGroup)
    byGroupSummary = t(apply(data,1,function(x){
        sapply(uniqueGroups,function(d) {
            subset = x[which(byGroup==d)]
                  do.call(myFunction,args=list(subset,...))
        })
    }))
    if(useNames) colnames(byGroupSummary) = uniqueGroups
    return(byGroupSummary)
}

# get list of unique (and non-equal) pairs within a vector

getUniquePairs <- function(X){
    
    Pairs = sapply(1:(length(X)-1), function(i) {
        a = sapply((i+1):length(X),function(j){
            d = X[c(i,j)]
        })
    })
    output <- matrix(unlist(Pairs), ncol = 2, byrow = T)
    return(output)
}


## function to order samples by number of occurances
order.by.number.occurrences <- function(group) {
    x = as.numeric(factor(group))
    x[is.na(x)] = max(x,na.rm=TRUE)+1
    order = sort(table(x)[x],decreasing=TRUE,index.return=TRUE)
    order = order$ix
    return(order)
}

################
# get powers of 10 format
get10power <- function(X,minPower=0,extra=NULL,...){
    Y = format(X,scientific=TRUE,...)
    parts=str_split(Y,"e")
    dec = sapply(parts,function(i) i[1])
    ten = sapply(parts,function(i) as.numeric(i[2]))
    if(!is.null(extra)) Ynew = sapply(1:length(Y),function(i) as.expression(bquote(.(extra)~.(dec[i])~x10^.(ten[i])))) else Ynew = sapply(1:length(Y),function(i) as.expression(bquote(.(dec[i])~x10^.(ten[i]))))
    Ynew[is.na(X)]=NA
    if(minPower>0) Ynew[which(X<(10^minPower))] = as.character(X[which(X<(10^minPower))])
    if(minPower<0) Ynew[which(X>(10^minPower))] = as.character(X[which(X>(10^minPower))])

    return(Ynew)
}

get10power2 <- function(X,minPower=0,digits=NULL){
    ten=log10(X)
    if(!is.null(digits)) ten = round(ten,digits)
    Ynew = sapply(1:length(X),function(i) as.expression(bquote(10^.(ten[i]))))
    Ynew[is.na(X)]=NA
    if(minPower>0) Ynew[which(X<(10^minPower))] = as.character(X[which(X<(10^minPower))])
    if(minPower<0) Ynew[which(X>(10^minPower))] = as.character(X[which(X>(10^minPower))])

    return(Ynew)
}

################
# Xtable helper functions

italic <- function(x){
paste0('{\\emph{', x, '}}')
}
large <- function(x){
paste0('{\\Large ', x, '}')
}
bold <- function(x){
paste0('{\\bfseries ', x, '}')
}

myAlign <- function(data,widths=NULL,totalWidth="\\textwidth",minCharWidth=3,nonNumeric=c(),rows=FALSE,returnWidths=FALSE,buffer=0.05){
# Get alignment settings for an xtable based on the input data. The same data must be used for the xtable() function!    
    nCols = ncol(data)
    # get max character width of columns.
    if(is.null(widths)) {
      if(!rows) widths = c(0,apply(data,2,function(x) max( max( nchar(as.character(x),type="width" )),minCharWidth)))
      if(rows) widths = c(apply(data,2,function(x) max( max( nchar(as.character(x),type="width" )),minCharWidth)))
    
    # any column that isn't numeric has smaller width contribution
      widths[1+nonNumeric] =  widths[1+nonNumeric]/3
      fractions = (1-buffer)*floor(100*widths/sum(widths))/100
      numericCols = apply(data,2,class)=="numeric"
    } else {
       if(!rows) fractions =c(0,widths)
       if(rows) fractions =c(widths)
    }
      alignment = paste0("|",paste(paste0(rep("p{",nCols+1),fractions,totalWidth,"}"),collapse="|"),"|")
    if(returnWidths) return(fractions) else return(alignment)
}


xdigits <- function (x, pad = TRUE, zap = getOption("digits")) {
    dig <- function(v) {
        if (is.numeric(v)) {
            v <- na.omit(v)
            v <- zapsmall(abs(v - floor(v)), zap)
            dec <- if (any(v > 0)) 
                max(nchar(v) - 2L)
            else 0L
        }
        else {
            dec <- 0L
        }
        return(dec)
    }
    is.2d <- length(dim(x)) == 2
    decimals <- if (is.2d) 
        sapply(as.data.frame(x), dig)
    else dig(x)
    output <- if (is.2d && pad) 
        c(0L, decimals)
    else decimals
    return(output)
}


#####################
# myHexbin.rect (like hexVP.abline but with rectangles)

myHexbin.rect <- function(hvp,x0,y0,x1,y1,
    col = "black", lty = 1, lwd = 2, ...){
  
    pushHexport(hvp, clip = "off")
    xx <- current.viewport()$xscale
    yy <- current.viewport()$yscale
    grid.lines(x=c(x0,x1),y = c(y0,y0), default.units = "native", gp = gpar(col = col.line, lty = lty, lwd = lwd))
    grid.lines(x=c(x1,x1),y = c(y0,y1), default.units = "native", gp = gpar(col = col.line, lty = lty, lwd = lwd))
    grid.lines(x=c(x1,x0),y = c(y1,y1), default.units = "native", gp = gpar(col = col.line, lty = lty, lwd = lwd))
    grid.lines(x=c(x0,x0),y = c(y1,y0), default.units = "native", gp = gpar(col = col.line, lty = lty, lwd = lwd))
    
    popViewport()
}               
  
myHexbin.xylab <- function(hvp,xlab,ylab,lineHeight=(-2),gp = gpar(fontsize = 16),...){
    pushHexport(hvp, clip = "off")
    grid.text(xlab, y = unit(lineHeight, "lines"), gp = gp,rot = 90,...)
    grid.text(ylab, x = unit(lineHeight, "lines"), gp = gp,...)
    popViewport()
}

#####################
# Arrange a grid plot

my.grid.arrange <- function (..., as.table = FALSE, clip = TRUE, main = NULL, sub = NULL, 
    left = NULL, legend = NULL, newpage = F) 
{
    grid.draw(arrangeGrob(..., as.table = as.table, clip = clip, 
        main = main, sub = sub, left = left, legend = legend))
}


################
# map aspect ratios

spainAspect = 1.34669
contextAspect = 1.315389
#ALL_context_map$aspect.ratio

################
# Make rectangle polygon from bounding box of a spatial data frame
boxPolygon <- function(box,id="box"){
    pol = rbind(box[,1],c(box[1,1],box[2,2]),box[,2],c(box[1,2],box[2,1]),box[,1])
    Polygons(list(Polygon(pol)),ID=id)
}

makeBboxPolygon <- function(mapDat,box=NULL){
  if(is.null(box)) box = bbox(mapDat)
  pol = rbind(box[,1],c(box[1,1],box[2,2]),box[,2],c(box[1,2],box[2,1]),box[,1])
  spPol = SpatialPolygons(list(Polygons(list(Polygon(pol)),ID="box")),proj4string = CRS(proj4string(mapDat)))
  dat = as.data.frame(matrix(1,1,1)); rownames(dat) = "box"; colnames(dat)="dummy"
  spPolDf = SpatialPolygonsDataFrame(spPol,data=dat)
  return(spPolDf)
}

# find the appropriate tiles (1degreex1degree) for a given point
pickTiles = function(tileCodes,centre,pad,type="raster"){
  
  if(type=="raster") exampleMap = raster(paste0(tileCodes[tileCodes!="NULL"][[1]],'.tif'))
  if(type=="shape") exampleMap = readOGR(dirname(paste0(tileCodes[tileCodes!="NULL"][[1]])),layer=basename(paste0(tileCodes[tileCodes!="NULL"][[1]])))
  if(length(pad)==1) pad= c(pad,pad)
    
  myBox = makeBboxPolygon(exampleMap,box=rbind(c(centre[1]-pad[1],centre[1]+pad[1]),c(centre[2]-pad[2],centre[2]+pad[2])))
  theseTiles = sapply(colnames(tileCodes),function(lat){
    lat=as.numeric(lat)
    p = sapply(rownames(tileCodes),function(lon){
      lon=as.numeric(lon)
      o = gIntersects( makeBboxPolygon(exampleMap,box=rbind(c(lon,lon+1),c(lat,lat+1))), 
                       myBox )
    })
  })
  out = tileCodes[theseTiles]
  if(length(out)==0) print("oh dear, there seem to be no tiles for you in the tileCodes matrix...")
  return(out)
}


################
# Make spatialpointslist from set of coordinates

makeSpPoints <- function(coords,mapDat=espMapTotal){
   points <- SpatialPoints(coords,proj4string=CRS(proj4string(mapDat)))
   return(points)
}

################
# Make spatial polygons that make an arrow!

# get points for arrow polygon based on a pair of coordinates
getArrowCoords <- function(coords,angle=30,width,frac1=2.5,frac2=2){
        
        r1 = atan( diff(coords)[2]/diff(coords)[1] )
        r2 = 2*pi*angle/360
        
        # rotate coordinates by this amount, around 0 to get a flat arrow
        newCoords = coords%*%cbind(c(cos(r1),sin(r1)),c(-sin(r1),cos(r1)))
        
        # width can't be wider than length of arrow
        lengthArrow = sqrt(diff(coords)[1]^2 + diff(coords)[2]^2)
        if(width>lengthArrow*tan(r2)/2){
          print("WARNING: width greater than possible! making width half the maximum size instead.")
          width = lengthArrow*tan(r2)/2
        }
        h=width/2
        direction = -sign(diff(newCoords[,1]))
        print(direction)
        
        baseCoord1 = newCoords[1,]+c(0,h)
        baseCoord2 = newCoords[1,]-c(0,h)
        
        dis = frac1*h/tan(r2)
        arrowCoord1 = newCoords[2,]+c(direction*dis,frac1*h)
        arrowCoord2 = newCoords[2,]+c(direction*dis,-frac1*h)
        
        r3 = r2/frac2
        dis2 = tan(pi/2 - r2 - r3)*(frac1*h-h)
        
        arrowCoord3 = newCoords[2,]+c(direction*(dis-dis2),h)
        arrowCoord4 = newCoords[2,]+c(direction*(dis-dis2),-h)
      
        allPoints = rbind(baseCoord1,baseCoord2,arrowCoord4,arrowCoord2,newCoords[2,],arrowCoord1,arrowCoord3,baseCoord1)   
        
        r4 = 2*pi - r1
        allPointsRotated = allPoints%*%cbind(c(cos(r4),sin(r4)),c(-sin(r4),cos(r4)))
        
        return(list(allPointsRotated,allPoints))
    }

# note only one colour can be used to fill these polygons    
getArrowSpPolygons <- function(Coords,widths=NULL,angle=30,width=1,mainMap=espMapTotal,...){
  #Coords is a list of 2x2 matrices with start and end points of arrows
  if(is.null(widths)) widths = rep(1,length(coords))

    arrowPolygons = sapply(1:length(Coords), function(i){
        out = getArrowCoords(Coords[[i]],angle=30,width=widths[i],frac1=2.5,frac2=2)
        Polygons(list(Polygon(out[[1]])),ID=i)
    })
    
    Arrows = SpatialPolygons(arrowPolygons, proj4string=CRS(proj4string(mainMap)))
    return(Arrows)
}



################
# Make lines polygon from two sets of coordinates

linesBetweenPoints <- function( begin.coord,end.coord,proj4str = proj4string(espMuniMap) ){
  
  l <- vector("list", nrow(begin.coord))
  colnames(begin.coord) = colnames(end.coord) = c("X","Y")
  for (i in seq_along(l)) {
      l[[i]] <- Lines(list(Line(rbind(begin.coord[i, ], end.coord[i,]))), as.character(i))
  }
  out = SpatialLines(l,proj4string = CRS(proj4str))
  
  dat = as.data.frame(matrix(1,length(out),1)); colnames(dat)="dummy"
  spLineDf = SpatialLinesDataFrame(out,data=dat)
  
  return(out)
}


################
# Create lines data frame for plotting a map scale under a given map 
# will be based on the same extent as the input for mapDat

createMapScalePlot <- function(mapDat,offsetFraction=20,...){
  # output is an spplot
  plotWidthLon = (extent(mapDat)[2]-extent(mapDat)[1])
  plotWidthKm = getWidthMapKM(mapDat)
  offsetY = (extent(mapDat)[4]-extent(mapDat)[3])/20
  textY = extent(mapDat)[3]-offsetY; textX = extent(mapDat)[1] + (extent(mapDat)[2]-extent(mapDat)[1])/2
  text1 = list("sp.text", c(textX,textY), paste0(plotWidthKm," Km"), which = 1,cex=2)
  xStart = extent(mapDat)[1]; xEnd = extent(mapDat)[2]
  scaleLines = linesBetweenPoints(begin.coord = rbind( c(xStart,textY),
                                          c(xEnd,textY),
                                          c(xStart,textY-offsetY/2),c(xEnd,textY-offsetY/2)),
                     end.coord = rbind( c(textX-plotWidthLon/20,textY),
                                        c(textX+plotWidthLon/20,textY),
                                        c(xStart,textY+offsetY/2),c(xEnd,textY+offsetY/2)),proj4str = proj4string(mapDat))
  myScale = spplot(scaleLines,colorkey=FALSE,sp.layout=list(text1),...)

  return(myScale)

}

getWidthMapKM <- function(mapDat){
  mapDatMeters = spTransform(makeBboxPolygon(mapDat),CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#NOTE: meters transformation requires a centering point. x_0 and y_0, as given here.
  plotWidthKm = round((extent(mapDatMeters)[2]-extent(mapDatMeters)[1])/1000,1)
  return(plotWidthKm)
}  # This isn't quite correct...

# THis works better. Based on gcDestination()
makeDistanceScaleBar <- function(myBox,distanceLon=50,dist.units = "km",scaleHeight=NULL,padding=1/20,scaleCols =add.alpha(c("white","black"),0.8),...){
  
  model = str_split(str_split(projection(myBox),"\\+ellps\\=")[[1]][2]," ")[[1]][1] # should be the same as the proj4string
  
  scalePos = bbox(myBox)[,1] + padding*diff(t(bbox(myBox)))  # bottom left relative to centre of plot with 10% padding
  if(is.null(scaleHeight)) scaleHeight = diff(bbox(myBox)[2,])/40
    
  # bottom-right coordinates of scale bar
  bl1hor1 <- gcDestination(lon = scalePos[1], lat = scalePos[2], bearing = 90, dist = distanceLon*2,dist.units = "km",model = model)
  bl1ver1 <- gcDestination(lon = scalePos[1], lat = scalePos[2], bearing = 0, dist = distanceLon/20,dist.units = "km",model = model)
  bl1ver2 <- gcDestination(lon = bl1hor1[1], lat = bl1hor1[2], bearing = 0, dist = distanceLon/20,dist.units = "km",model = model)
  halfway1 = colMeans(rbind(bl1hor1,scalePos)) # half way between points
  halfway2 = colMeans(rbind(bl1ver1,bl1ver2)) # half way between points (note doesn't take curvature into account)
    
  # make polygons
  pol1 = rbind(scalePos,halfway1,halfway2,bl1ver1,scalePos)
  pol2 = rbind(halfway1,bl1hor1,bl1ver2,halfway2,halfway1)
  
  spPol = SpatialPolygons(list(Polygons(list(Polygon(pol1)),ID="box1"),
                               Polygons(list(Polygon(pol2)),ID="box2")),
                          proj4string = CRS(proj4string(myBox)))
  #dat = as.data.frame(matrix(1,2,1)); rownames(dat) = c("box1","box2"); colnames(dat)="dummy"
  #spPolDf = SpatialPolygonsDataFrame(spPol,data=dat)
  spPols =  list("sp.polygons",spPol,col="black",fill=scaleCols,which=1,lwd=2)
  
  # add text
  text1 = list("sp.text",bl1ver1+c(0,scaleHeight/1.7),"0",which=1,...)
  text2 = list("sp.text",bl1ver2+c(0,scaleHeight/1.7),paste0(distanceLon*2," ",dist.units),which=1,...)
  
  return(list(spPols,text1,text2))
}

makeDistanceScaleLine <- function(myBox,distanceLon=50,dist.units = "km",bearing=90,padding=1/20,scaleCols =add.alpha(c("white","black"),0.8),...){
  
  model = str_split(str_split(projection(myBox),"\\+ellps\\=")[[1]][2]," ")[[1]][1] # should be the same as the proj4string
  
  scalePos = bbox(myBox)[,1] + padding*diff(t(bbox(myBox)))  # Starting position is bottom left relative to centre of plot with 10% padding
  if(is.null(scaleHeight)) scaleHeight = diff(bbox(myBox)[2,])/40
    
  # right-hand coordinates of line
  bl1hor1 <- gcDestination(lon = scalePos[1], lat = scalePos[2], bearing = bearing, dist = distanceLon*2,dist.units = "km",model = model)
    
  # make lines
  Line = linesBetweenPoints(begin.coord = scalePos,end.coord = bl1hor1)
  Line = list("sp.lines",Line,...)
  
  # add text
  text1 = list("sp.text",scalePos+c(0,scaleHeight/2),"0",which=1,...)
  text2 = list("sp.text",bl1hor1+c(0,scaleHeight/2),paste0(distanceLon*2," ",dist.units),which=1,...)
  
  return(list(Line,text1,text2))
}

################
# FUnction to jitter points (also a copy of this is used in /samples/SPAIN.A.geocode.R)

Jitter <- function(x,y,jitterSize,myMap=espMapTotal,above=FALSE){
  if ((is.na(x)==T) | (x=='character(0)')) {return(c(0,0))} else {
   test <- over(SpatialPoints(t(c(x,y)),proj4string=CRS(proj4string(myMap))),myMap)
   if(is.na(test[1])){ return(c(x,y)) } else {
   test <- NA
  while (is.na(test[1])==T) {
  xhat <- jitter(0,amount=jitterSize)
  xpos <- x + xhat
  if(!above) ypos <- y + runif(1,-sqrt(jitterSize^2-xhat^2),sqrt(jitterSize^2-xhat^2))
  if(above) ypos <- y + runif(1,0,sqrt(jitterSize^2-xhat^2)) # only jitter to semi-circle above original point

  #ypos <- jitter(y,amount=jitterSize)
  test <- over(SpatialPoints(t(c(xpos,ypos)),proj4string=CRS(proj4string(myMap))),myMap)         
   }
   return(c(xpos,ypos))
  }
 }
 }
 

################
# FUnction to jitter points nicely within a particular mapping window

getNicePointsLayout <- function(indSubset,max.delta=0.0001,myBox,rad=NULL){
     #max.delta sets maximum distance that a node can move from it's original position
     # indSubset is a dataFrame derived from stuff$D
  
      indSubset$Xmuni.ave.grand.custom2=indSubset$Xmuni.ave.grand.custom
      indSubset$Ymuni.ave.grand.custom2=indSubset$Xmuni.ave.grand.custom
      
      # exclude anyone whose original geocode is not within the bounds of myBox
      exc = ( indSubset$Xmuni.ave.grand.precise < bbox(myBox)[1,1] ) | ( indSubset$Xmuni.ave.grand.precise > bbox(myBox)[1,2] ) | ( indSubset$Ymuni.ave.grand.precise < bbox(myBox)[2,1] ) | ( indSubset$Ymuni.ave.grand.precise > bbox(myBox)[2,2] )
      print("exclusions outside box:")
        print(sum(exc))
        
      indSubset = indSubset[!exc,]
      
      fixedNodes = unique(indSubset[,c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")])
        
      EDGES = NODES = GRAPHS = list()
      maxDeltaPoints = c()
      
      for(i in 1:nrow(fixedNodes)){
        print(i)
        x = fixedNodes[i,]
        rownames(x) = paste0("node",i)
        inds = which((indSubset[,c("Xmuni.ave.grand.precise")]==x[1,1]) & (indSubset[,c("Ymuni.ave.grand.precise")]==x[1,2]))
            jitteredNodes = indSubset[inds,c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")] + rmvnorm(length(inds),mean = c(0,0),sigma = diag(c(0.0001,0.0001))) # randomly perturb a tiny bit
            Nodes = as.data.frame(t(cbind(t(x),t(jitteredNodes))),stringsAsFactors=FALSE); Nodes$name = rownames(Nodes)
            Nodes = Nodes[,c(3,1,2)]
            Edges = cbind(rep(paste0("node",i),dim(jitteredNodes)[1]),rownames(jitteredNodes))
            g = graph.data.frame(d = Edges ,directed = FALSE,vertices = Nodes)
            EDGES[[i]] = Edges
            NODES[[i]] = Nodes
            GRAPHS[[i]] = g
        
            # set max delta proportional to the number of inds in the group ()
            maxDeltaPoints = c( maxDeltaPoints,rep(max.delta*(1 + log10(length(inds))),dim(Nodes)[1]) )

        if(length(inds)>2){
          print(i)
        
            l = layout.lgl(g,root=paste0("node",i));
            newPositions = l; rownames(newPositions) = V(g)$name
            newPositions[,1] = newPositions[,1] -  l[1,1] + fixedNodes[i,1];
            newPositions[,2] = newPositions[,2] -  l[1,2] + fixedNodes[i,2];
            newPositions = newPositions[-1,]
            indSubset[rownames(newPositions),"Xmuni.ave.grand.custom2"] = newPositions[,1]
            indSubset[rownames(newPositions),"Ymuni.ave.grand.custom2"] = newPositions[,2]
        }
      }
      
      myNodes = rbind_all(NODES); myEdges = abind(EDGES,along=1)
      
      qEdges = qgraph(myEdges,DoNotPlot=TRUE); conc=unique(cbind(c(qEdges$Edgelist$from,qEdges$Edgelist$to),c(myEdges[,1],myEdges[,2])))
      NodeCoords = as.matrix(myNodes[match(conc[,2],myNodes$name),2:3]); 
      fixedNodes = grepl("node",conc[,2])
      myFixedNodes = NodeCoords; myFixedNodes[!fixedNodes,]=c(NA,NA)
      colors = rep("black",length(conc)); colors[fixedNodes]="red"
      layoutNice = qgraph.layout.fruchtermanreingold(cbind(qEdges$Edgelist$from,qEdges$Edgelist$to),repulse.rad = rad,
                                                     init = NodeCoords,max.delta=maxDeltaPoints,
                                                     vcount=dim(NodeCoords)[1],constraints=myFixedNodes)
      q1 = qgraph(myEdges,layout=layoutNice,layout.orig=layoutNice,DoNotPlot=TRUE,rescale=TRUE,color=colors,vsize=1)
      #plot(q1)
      #q2 = qgraph(myEdges,layout=NodeCoords,layout.orig=NodeCoords,DoNotPlot=TRUE,rescale=TRUE,color=colors,vsize=1)
      #plot(q2)

      print("getting data frames etc.")
      # plot with boxplot
      rownames(layoutNice) = conc[,2];
      test = makeSpPoints(coords = layoutNice[!fixedNodes,],mapDat=myBox)
      newJitterLines = linesBetweenPoints(begin.coord=indSubset[rownames(layoutNice[!fixedNodes,]),c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")],
                         end.coord=layoutNice[!fixedNodes,])
      
      origPoints = list("sp.points",makeSpPoints(indSubset[,c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")]),
                  cex=1,col="lightgray",pch=16)

      indSubset[rownames(layoutNice[!fixedNodes,]),c("Xmuni.ave.grand.custom2","Ymuni.ave.grand.custom2")] = layoutNice[!fixedNodes,]
      
      return(list( layoutNice,newJitterLines,origPoints,indSubset ))
}


################
# flit matrices for images

imgflip<-function(x) {t(x[nrow(x):1,])}

################
# plot themes

# lattice plots
theme.novpadding <-
   list(layout.heights =
        list(top.padding = 0,
       main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	   xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 3,
 	    right.padding = 0
        ))
             

#levelplot(raster(sum),col.regions=some.colorsEnd2,colorkey=T,scales=list(draw = FALSE),margin=F,
#                              par.settings = theme.novpadding)

################
# shapes

sh <- c(25,21,22,23,24)
numToShape <- rbind(seq(0,(length(sh)-1)),sh)
colnames(numToShape) <- numToShape[1,]


################
# colours

getNumToColour <- function(split.matrix,maxSplit=ncol(split.matrix),valueShift=0){
  tab <- table(split.matrix[,maxSplit])
  splitOrder <- unique(split.matrix[,maxSplit])
  tabSort <- tab[as.character(splitOrder)]
  x <- cumsum(tabSort)/sum(tabSort)
  y <- x
  y[which(x>valueShift)] <- x[which(x>valueShift)]-valueShift
  y[which(x<valueShift)] <- 1-(valueShift-x[which(x<valueShift)])
  co <- hsv(y,rep(c(0.6,0.8,1,1,1),length(tab)/5),rep(c(1,1,0.8,0.6,0.4),length(tab)/5))
  splitCols <- rbind(as.character(splitOrder),co)
  splitCols <- splitCols[,order(splitOrder)]
  colnames(splitCols) <- splitCols[1,]
 pizza(splitCols[2,])
  return(splitCols)
}

distinctColours <- read.table('~/Documents/ClareDPhil/DPhil/Spain/Other/distinctColours-Intense-60.txt',stringsAsFactors=F, header=F)[,1]

getNumToColour2 <- function(split.matrix,maxSplit=ncol(split.matrix),colourList=distinctColours,random=F){
  
  tab <- table(split.matrix[,maxSplit])
  splitOrder <- unique(split.matrix[,maxSplit])
  tabSort <- tab[as.character(splitOrder)]
  if(length(splitOrder)<=length(colourList)){
  if(random==T) co <- sample(c(1:7,colourList),length(splitOrder),replace=F) else co <- c(2:7,colourList)[1:maxSplit]
  splitCols <- rbind(as.character(splitOrder),co)
  splitCols <- splitCols[,order(splitOrder)]
  colnames(splitCols) <- splitCols[1,]
    pizza(splitCols[2,])
    return(splitCols)
  } else {print('cant make colours as need more colours in colourList')}
}

getNumToColourDistant <- function(split.matrix,maxSplit=ncol(split.matrix),valueShift=0){
  tab <- table(split.matrix[,maxSplit])
  splitOrder <- unique(split.matrix[,maxSplit])
  tabSort <- tab[as.character(splitOrder)]
  
  coDist <- c()
  index <- 1
  c3 <- "red"
  for (i in splitOrder){
    c <- splitComp(c3,plot=F)
    c2 <- c[2]
    c3a <- rgb2hsv(col2rgb(c[3]))
    h <- c3a[1,1]
    s <- c3a[2,1]
    v <- c3a[3,1]
    #c3a[,1] <- c3a + c(runif(1,h+0.03,1-h),runif(1,0.7-s,1-s),runif(1,0.3-v,1-v))
    c3a[,1] <- c3a + c(0.03,runif(1,0.8-s,1-s),runif(1,0.3-v,1-v))
    
    if(c3a[1,1]>1){
       c3a[1,1] <- c3a[1,1]-1
    }
    #print(c3a)
    c3 <- hsv(c3a[1,1],c3a[2,1],c3a[3,1])
    if(index==1) {
      c2 <- "red"
      c3 <- "red"
      index <- 2
    }
    coDist <- c(coDist,c2)
    
  }
  
  splitCols <- rbind(as.character(splitOrder),coDist)
  splitCols <- splitCols[,order(splitOrder)]
  colnames(splitCols) <- splitCols[1,]
 pizza(splitCols[2,])
return(splitCols)
}

# Function to find a set of n colours maximally distant in the 3-dimesional hsv colour space, based on a given vector of colours. Default is colors()
getColoursDistant <- function(n,colorSet=colors()){
    # matrix where each color is a row of three elements
    colorsHsv=t(rgb2hsv(col2rgb(colorSet)))
    if(dim(unique(colorsHsv))[1]<dim(colorsHsv)[1]) warning(paste0("Your colorSet has non-unique hsv values."))
    colorsHsv = unique(colorsHsv)
    nCols = dim(colorsHsv)[1]
    print(paste0(nCols," unique colors."))
    
    if( nCols <= n ) {
      warning(paste0("Your colorSet has ",n-nCols," too few colors. Returning the whole set."))
      return(colorSet)
    }
    
    # compute matrix of euclidian distances and find a maximally distant set
    subset <- colorsHsv
    alldist <- as.matrix(dist(subset))
    SVD = svd(alldist)
    
    h = sapply(1:n,function(x) {
      y = SVD$u[,x]
      c(which(y==min(y)),which(y==max(y)))
    })
    sub = c()
    i=1
    while(length(sub)<n){
      sub = c(sub,h[,i][!h[,i]%in%sub])
      i=i+1
    }
    subset = colorsHsv[sub,]
    
   # while (nrow(subset) > n) {
  #      #cdists = rowMeans(alldist)
  #      diag(alldist)=NA
  #      close <- which(cdists == min(cdists))
        # remove the one with the next smallest distance
  #      closest = close[order(apply(alldist[close,-close],1,min,na.rm=TRUE))][1]
  ####      cdists = apply(alldist,1,min,na.rm=TRUE)
  #      close = c(close,closest)
  #      subset <- subset[-closest,]
  #      alldist <- alldist[-closest,-closest]
  #  }
    # convert back to rgb
    out = apply(subset,1,function(x) hsv(x[1],x[2],x[3]))
    #pizza(out)
    return(out)
}

# try something else: find maximum possible minimum distance between a set of n points.
isFeasibleDOESNTWORK <- function(dist, k, distMat,dontFindMin=TRUE){
    # Can we construct a subset at least as big as k such that no value is less than dist?  
    #   newMat = rbind(c(0,0,0,0,1,0),c(0,0,1,0,0,1),c(0,0,0,0,0,0),c(1,0,0,0,0,0),c(0,0,0,0,0,1),c(0,0,1,0,0,0))
    #  newMat[newMat==0]=2; newMat[newMat==1]=0;newMat[newMat==2]=1
 
    newMat = distMat; newMat[distMat<dist]=0; newMat[distMat>=dist]=1
    diag(newMat)=1
    
    toKeep = which(rowSums(newMat)>k)
    subMat = newMat[toKeep,toKeep]
    image(subMat)
    hMat = subMat;colnames(hMat)=1:ncol(hMat); hMat[]=0; wMat=hMat; area_max=0
    
    for ( x in 1:ncol(hMat) ){
      for ( y in 1:ncol(hMat) ){
        if(subMat[x,y]==0) next
        if(x==1) hMat[x,y] = 1 else hMat[x,y] = 1 + hMat[x-1,y]
        if(y==1) wMat[x,y] = 1 else wMat[x,y] = 1 + wMat[x,y-1]
        minw = wMat[x,y]
        
        for( dh in 1:hMat[x,y] ){
          if(x!=y) next # enforce symmetry on rows and columns
          if(dh < k) next
          minw = min(minw, wMat[x-dh+1,y])
          if(dh != minw) next
          area = dh*minw
          print(c(x-1,y-1,dh,minw))
          if ( area > area_max[1] ){
                area_max = c(area, x-dh+1, y-minw+1, x, y)
                if(dontFindMin) break
          }
          
        }
      }
    }
    
    if( sum(subMat[area_max[2]:area_max[4],area_max[3]:area_max[5]]!=1)!=0 ) print("SOMETING WONG!")
    distMat[rownames(outMat),colnames(outMat)]
    
    
}    
    
isFeasible <- function(dist, k, distMat,returnRows=FALSE){
    # Can we construct a subset at least as big as k such that no value is less than dist?  
    #   newMat = rbind(c(0,0,0,0,1,0),c(0,0,1,0,0,1),c(0,0,0,0,0,0),c(1,0,0,0,0,0),c(0,0,0,0,0,1),c(0,0,1,0,0,0))
    #  newMat[newMat==0]=2; newMat[newMat==1]=0;newMat[newMat==2]=1
 
    newMat = distMat; newMat[distMat<dist]=0; newMat[distMat>=dist]=1
    diag(newMat)=1
    toKeep = which(rowSums(newMat)>k)
    if(length(toKeep)<k) return(FALSE)
    subMat = newMat[toKeep,toKeep]
    
    a = subMat
    for(i in 1:ncol(subMat)){
      a = a[order(a[i,]),order(a[i,])]
    }
    #image(a)
    #image(subMat)
    # find maximum contiguous block
    for( i in  1:(nrow(a)-k)){
      allOnes = sum(a[i:(i+k-1),i:(i+k-1)]==0)==0 
      if(allOnes) {
        if(returnRows){
          out = colnames(a)[i:(i+k-1)]
          return(out)
        } else {
          return(allOnes)
        }
      }
    }
    return(FALSE)
}
solveMaxMin <- function(k, distMat,fixed=NULL,tol=NULL){
    low = min(distMat,na.rm=TRUE) #// definitely small enough.
    high = max(distMat,na.rm=TRUE) + 1 #// definitely too big
    if(is.null(tol)) tol = (high-low)/1000
    while ( (high - low) > tol ){
        mid = (low + high) / 2
        testFeasible = isFeasible(mid, k, distMat)
        if (testFeasible){
            low = mid
            print("yay, keep going higher!")
        } else {
            print("Nope, switch back to something smaller")
            high = mid
        }
    }
   getFeasible = isFeasible(low, k, distMat,returnRows=TRUE)
   print(low)
   return(as.numeric(getFeasible))
} 

getColoursDistant2 <- function(cols=colors(),n=10){
  
    colorsHsv=t(rgb2hsv(col2rgb(cols)))
    colorsHsv = unique(colorsHsv)
    colorsHsvSort = colorsHsv[order(colorsHsv[,1]),]
    colorsHsvSort = colorsHsvSort[order(colorsHsvSort[,2]),]
    colorsHsvSort = colorsHsvSort[order(colorsHsvSort[,3]),]
    
    distMat = as.matrix(dist(colorsHsvSort))
    diag(distMat)=NA
    subset = solveMaxMin(n,distMat)
    o = apply(colorsHsv[as.numeric(subset),],1,function(x) hsv(x[1],x[2],x[3]))
    return(o)
    #pizza(o)
}

# This one is actually good!
getColoursDistant3 <- function(n,huFrac=0.6,startHue=0,minLum=20,maxLum=80,chroma=10){
  if(huFrac<1/n) huFrac=1/n
  h=seq(360/n+startHue,360+startHue,by=360/n)
  h=h%%360
  l=seq((maxLum-minLum)/(n*huFrac),maxLum,by=(maxLum-minLum)/(n*huFrac))
  print(l)
  myCols = sapply(1:n,function(x) {
        hcl(h=h[x],l=l[x%%length(l)+1],c=chroma,fixup=TRUE)
      })
  return(myCols)
  #pizza(myCols)
}

################
#latin1 encoding concordance with ASCII
ASCII.Latin1Chart <- read.delim('~/Documents/ClareDPhil/DPhil/Spain/Other/Latin1-ASCII-chart.txt',header=F,fileEncoding="Latin1",stringsAsFactors=F)
rownames(ASCII.Latin1Chart) <- ASCII.Latin1Chart$V2

################
# stuff for chromopainter output
#matrixTypes<- c("dataraw","chunklengthsmatrix","meanchunklengthsmatrix")
matrixTypes<- c("chunkcountsmatrix","chunklengthsmatrix","meanchunklengthsmatrix")
names(matrixTypes) <- c("ChunkCounts","ChunkLengths","meanChunkLengths")

################
# convert colour to kml
col2kml <- function(colorString){
  color=strsplit(col2hex(colorString),"")[[1]][-1]
  paste("ff",
       paste(color[5:6],collapse=""),
       paste(color[3:4],collapse=""),
       paste(color[1:2],collapse=""),
                   sep="")
}

#####################
# Plotting PCA of SPANISH data
#####################
  
 myPCA <- function(data) {
   
   D <- apply(data,MARGIN=2,function(x) (x-mean(x))/sd(x))
SVD <- svd(t(D))
SVD.PCs <- as.data.frame(t(t(SVD$v)*SVD$d)) # scale PCs by eigen values
SVD.vals <- as.data.frame(SVD$d)
SVD.loads <- as.data.frame(SVD$u)
rownames(SVD.PCs) <- rownames(D)
  names(SVD.PCs) <- paste("PC",c(1:ncol(SVD.PCs)),sep="")
   return(list("SVD.PCs"=SVD.PCs,"SVD.vals"=SVD.vals,"SVD.loads"=SVD.loads))
 }


# basic PCA function

basicPCAplot <- function(PCdata,PCX,PCY,colVar,legend=NULL,eig_val=NULL,pointSize=2,alpha=1,title="",legend_title=NULL,borders=F,legCols=1,...){

  .e <- environment()
  x <- PCdata
  x$geog <- PCdata[,colVar]
    
 if(is.null(legend)){
   legend <- makeFactorLegend(x$geog)
 }
  
  colors <- legend$colour
  names(colors) <- rownames(legend)
  borderColors <- legend$colour
  if("border"%in%colnames(legend))  borderColors <- legend$border else borderColors <- legend$colour
  names(borderColors) <- rownames(legend)
  if (borders==T) borderColors <- rep("black",nrow(legend))
  shapes <- as.numeric(legend$shape)
  names(shapes) <- rownames(legend)
  labels <- legend$Factor
  print(labels)
  names(labels) <- rownames(legend)
  if(is.null(legend_title)){
  g <- guide_legend(colVar,ncol=legCols)
  } else {
    g <- guide_legend(legend_title,ncol=legCols)
  }
  
  if (is.null(eig_val)==FALSE) {
  PCXval <- round(100*eig_val[as.numeric(str_extract_all(PCX,"[[:digit:]]{1,}"))]/sum(eig_val),1)
  PCYval <- round(100*eig_val[as.numeric(str_extract_all(PCY,"[[:digit:]]{1,}"))]/sum(eig_val),1)
  xName <- paste(PCX," (",PCXval,"%)",sep="")
  yName <- paste(PCY," (",PCYval,"%)",sep="")
  } else {
    xName <- PCX
    yName <- PCY
  }
  print(c(PCX,PCY))
  
  
  plot <- ggplot(x,  aes(x=x[,PCX], y=x[,PCY]),environment=.e) + 
    #geom_point(aes(col=geog,shape=geog,fill=geog),size=pointSize,alpha=alpha,size=pointSize) +
    geom_point(aes(col=geog,shape=geog,fill=geog),size=pointSize,size=pointSize) +
    scale_color_manual(values=borderColors,labels=labels) +
    scale_shape_manual(values=shapes,labels=labels) +
    scale_fill_manual(values=colors,labels=labels) +
    guides(colour=g,shape=g,size=g,fill=g) +
    ggtitle(title) +
    scale_x_continuous(name=xName) +
    scale_y_continuous(name=yName) +
    theme(plot.title=element_text(size=rel(2)),legend.text=element_text(size=rel(0.8)),legend.title=element_text(size=rel(1.2)),plot.background=element_rect(fill="transparent"),...)
  
  return(plot)
}



SPAIN.PCA.plot <- function(PC,set="A",col_var="geog",title="PCA plot",PCX="PC1",PCY="PC2",ByColours="geography",
                           split=NULL,other_name=NULL,eig_val=NULL,alpha=0.9,pointSize=3,...) {

  # set up plot datasets (subset legend.master to get correct colourkey)
  g <- guide_legend("Espana comunidad")
  if(!is.null(col_var)) PC$geog <- PC[,col_var]
  
  if ((set=="A")&&(ByColours=="geography"))  {
        #PC$geog <- unlist(SPAIN.A.attributes.all[,"geog2"])
  legend <- legend.master[(legend.master$Freq.A>0)&(!is.na(legend.master$Freq.A)),]
  legend$labels = legend$SPAIN.A.name
  rownames(legend)<- legend$SPAIN.A.name
   g <- guide_legend("Autonomous Community")
  }
  
  if ((set=="Context")&(ByColours!="other")){
    c = ByColours
  legend <- legend.master.context[which(legend.master.context[,c]!="Spain"),]
  legend <- rbind(legend,legend.master.context[1,])
  legend$labels = legend[,c]
    rownames(legend)<- legend[,c]
    g <- guide_legend(c)
    
   }
  
  
   if (ByColours=="hospitals")  {
  
  legend <- legend.hospitals
  g <- guide_legend("hospital")
        
  }
  
  
    if (ByColours=="other")  {
        #PC$geog <- unlist(SPAIN.A.attributes.all[,"geog2"])
  
  legend <- as.data.frame(unique(PC$geog),stringsAsFactors=FALSE)
  colnames(legend) <- "labels"
  if (nrow(legend)>2) legend$colors <- c(hsv(seq(0,1,length.out=(nrow(legend)-1)),s=c(0.7,1),v=c(0.7,0.8)),"black")
  if (nrow(legend)<=2) legend$colors <- hsv(seq(0,1,length.out=nrow(legend)))
  legend$shapes <- c(1:nrow(legend))%%6
  if (length(which(is.na(legend$labels)))>0) legend$labels[which(is.na(legend$labels))] <- "na"
  rownames(legend)<- legend$labels
  g <- guide_legend(paste(col_var))
        
  }
  
  if ((set=="A&C")&&(ByColours=="geography"))  {
        #PC$geog <- unlist(SPAIN.A.attributes.all[,"geog2"])
  legend1 <- legend.master[(legend.master$Freq.A>0)&(!is.na(legend.master$Freq.A)),]
  legend1$labels = legend1$SPAIN.A.totals
  rownames(legend1)<- legend1$SPAIN.A.name
   legend.master.context$shapes=1
  names(legend.master.context)[1:2] <- c("labels","colors")
  legend <- rbind(legend1[,c("labels","colors","shapes")],legend.master.context)
  g <- guide_legend("Autonomous Community")
  }
  
  if ((set=="B")&&(ByColours=="geography")) {
   #   if ("geog"%in%names(PC)){ # do nothing
  #  } else {
    #PC$geog <- SPAIN.B.attributes.all$Comunidad[SPAIN.B.attributes.all$excluded==0]
   # }
    legend <- legend.master[(legend.master$Freq.B2>0)&(!is.na(legend.master$Freq.B)),]
    legend$labels = legend$SPAIN.B.totals
    rownames(legend) <- legend$SPAIN.B.name
    g <- guide_legend("Autonomous Community")
  }
  
  if (ByColours=="split"){
   PC$geog <- PC[,split]
  legend <- numToColour[2,unique(PC[,split])]
  legend <- as.data.frame(legend)
  legend$labels <- as.character(unique(PC[,split]))
  names(legend)[1] <- "colors"
   legend$colors <- as.character(legend$colors)
  legend$shapes <- as.numeric(legend$labels)%%6
  rownames(legend) <- legend$labels
  
  g <- guide_legend("fineSTRUCTURE cluster")
  PC$geog <- factor(PC$geog,levels=unique(PC$geog))
          
      }
  
  # Get eigen value info and calculate '% variance explained' for each PC axis
  if (is.null(eig_val)==FALSE) {
  PCXval <- round(100*eig_val[as.numeric(str_extract_all(PCX,"[[:digit:]]{1,}"))]/sum(eig_val),1)
  PCYval <- round(100*eig_val[as.numeric(str_extract_all(PCY,"[[:digit:]]{1,}"))]/sum(eig_val),1)
  xName <- paste(PCX," (",PCXval,"%)",sep="")
  yName <- paste(PCY," (",PCYval,"%)",sep="")
  } else {
    xName <- PCX
    yName <- PCY
  }
  
  .e <- environment()
  x <- PC
  
  
  colors <- legend$colors
  names(colors) <- rownames(legend)
  shapes <- legend$shapes
  names(shapes) <- rownames(legend)
  labels <- legend$labels
  names(labels) <- rownames(legend)
  
  plot <- ggplot(x, aes(x=x[,PCX], y=x[,PCY]),environment=.e) + geom_point(aes(col=geog,shape=geog),alpha=alpha,size=pointSize) +
    scale_color_manual(values=colors,labels=labels) +
    scale_shape_manual(values=shapes,labels=labels) +
    guides(colour=g,shape=g,size=g) +
    ggtitle(title) +
    scale_x_continuous(name=xName) +
    scale_y_continuous(name=yName) +
    theme(plot.title=element_text(size=rel(2)),legend.text=element_text(size=rel(1)),legend.title=element_text(size=rel(1.2)))
  
  return(plot)
}


#Function for plotting PCs (SPAIN A for now) 
####### Input: PCA data of the type used in SPAIN.PCA.plot

plotPCsbyHospitalandGeog <- function(inputdata,eig_val,filename,num_pcs=q,set="A",plotType="grid"){
  
if(plotType=="grid"){
  pdf(filename,width=12,height=10)
     D=3
     H=3
 layoutMat <- matrix(c(1:nrow(legend.hospitals)),nrow=D,ncol=H)
   
  for (j in c(2:num_pcs)) {
 
  grid.newpage()
     inner=viewport(0.501,0.501,height=0.95,width=0.95,layout = grid.layout(D,H))
     outer=viewport(width=1,height=1)
    pushViewport(outer)
  pc1 <- paste("PC",1,sep="")
  pc2 <- paste("PC",j,sep="")
  PCXval <- round(100*eig_val[as.numeric(str_extract_all(pc1,"[[:digit:]]{1,}"))]/sum(eig_val),1)
  PCYval <- round(100*eig_val[as.numeric(str_extract_all(pc2,"[[:digit:]]{1,}"))]/sum(eig_val),1)
      
      # x,y labels
      print(grid.text(paste(pc1," (",PCXval,"%)",sep=""),x=0.5,y=0.02))
      print(grid.text(paste(pc2," (",PCYval,"%)",sep=""),x=0.015,y=0.5,rot=90))
      
  
  
     pushViewport(inner)
 
     for(K in c(1:nrow(legend.hospitals))){
        k <- unique(inputdata$hospital)[K]
        Row=which(layoutMat==K,arr.ind=T)[1]
        Col=which(layoutMat==K,arr.ind=T)[2]

   inputdata$hospital2 <- "black"
  inputdata$hospital2[inputdata$hospital!=k] <- NA
  
  
   p1<- SPAIN.PCA.plot(PC=inputdata,set=set,PCX=pc1,PCY=pc2,split=split,ByColours="geography",
                        title=waiver(),alpha=0.3,eig_val=eig_val,pointSize=2)
   
   p <- p1 + theme(panel.margin=unit(rep(0.01,4),"npc"),
                   plot.margin=unit(rep(0,4),"npc"),
                   axis.text=element_text(size=10),
                   axis.title=element_blank(),
                   text=element_text(size=5),
                   legend.position = "none",
                    plot.background=element_rect(fill="transparent"),
                    panel.border=element_blank()) +
         geom_point(color=inputdata$hospital2,shape=4,cex=1) + ggtitle(k)
        
    if ((Col==1)&(Row!=D)) p <- p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    if ((Col!=1)&(Row!=D)) p <- p + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
    if ((Col!=1)&(Row==D)) p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
        
  print(p, vp = viewport(layout.pos.row = Row, layout.pos.col = Col))
    
          }

        }
  
    dev.off()

}
  
if (plotType=="single"){
  for (j in c(2:num_pcs)){
      pc1 <- paste("PC",1,sep="")
      pc2 <- paste("PC",j,sep="")  
      pdf(paste(sub(".pdf","",filename),pc2,"pdf",sep="."),width=12,height=10)
     
          for(K in c(1:nrow(legend.hospitals))){
                        k <- unique(inputdata$hospital)[K]
                  inputdata$hospital2 <- "black"
                  inputdata$hospital2[inputdata$hospital!=k] <- NA
                        
                  p1<- SPAIN.PCA.plot(PC=inputdata,set=set,PCX=pc1,PCY=pc2,split=split,ByColours="geography",
                                        title=waiver(),alpha=0.3,eig_val=eig_val)
                  p <- p1 + geom_point(color=inputdata$hospital2,shape=4,cex=3) + ggtitle(k)
                  print(p) 
            }
      dev.off()
  }
}
 
 }



#####################
# Plotting fineSTRUCTURE output
#####################

#Function for fixing row and column names - removes leading 'X' put on by R 
####### Input: data of type 'chunkcounts'

    fixNames <- function(x) {
        rownames(x) <- sub("X","",rownames(x))
        colnames(x) <- sub("X","",colnames(x))
        return(x)
      }

#Function for making a colour legend for arbitrary factors
####### Input: a vector made up of some factors
makeFactorLegend <- function(f,numeric=F,newColours=NULL,shapes=NULL){
    S <- unique(f)
    S <- S[!is.na(S)]
    t <- as.numeric(as.factor(S))
    if (is.null(newColours)) c <- numToColour[2,t] else c <- newColours[1:length(t)]
    if(is.null(shapes)) s <- numToShape[2,as.character(t%%5)]
    if(!is.null(shapes)){
      s <- shapes[1:length(t)]
    }
      
  context_cols <- as.data.frame(cbind(S,c,as.vector(s)),stringsAsFactors=F)
  rownames(context_cols) <- as.character(context_cols[,1])
  colnames(context_cols) <- c("Factor","colour","shape")
  context_cols$shape <- as.numeric(context_cols$shape)
  return(context_cols)
}


#Function for making nice labels when many labels are the same
####### Input: data of type - vector of labels where length(unique(<vector>)) < length(<vector>). 
# The vector should also be sorted so that the same values appear in the same block.

niceLabels <- function(x) {
        if (length(unique(x)) == length(x)) {
            print("no need to do this as there won't be any change")
        }
        
        labels <- rle(x)
        if (length(labels$values)>length(unique(x))) print("NOTE: duplicated labels aren't contiguous in labels vector")

        nice <- unlist(sapply(c(1:length(labels$values)),function(i){
                n = labels$lengths[i]
                y = labels$values[i]
                pos = round(n/2)
                out = c(rep("",pos),y,rep("",n-pos-1))
                #nice <- c(nice,out)
                return(out)
          }))
        names(nice) <- NULL
        return(nice)
}

#Function for plotting convergence of fineSTRUCTURE MCMC
####### Input: data of type 'mcmcdata'
 
plotMCMCParameters <- function(data,title) {
  
a <- ggplot(data,aes(x=Number,y=Posterior)) + geom_point() + 
        scale_y_continuous(labels = comma) + ggtitle(title) + 
      xlab("thinned MCMC sample") + ylab("Log-posterior") + theme(plot.background=element_rect(fill="transparent"))
b <- ggplot(data,aes(x=Number)) +      
        geom_point(aes(y=AccHyper,colour="blue"))+ 
        geom_point(aes(y=AccIndiv,colour="red"))+ 
        geom_point(aes(y=AccSAMS,colour="green")) +
        geom_point(aes(y=AccMS,colour="purple")) +
        scale_color_manual(values=c("blue","red","green","purple"),
                           labels=c("Hyper","Indiv","SAMS","MS"),name=NULL)+
        ylab("Acceptance probability") + xlab("thinned MCMC sample")+
        ggtitle(title) + theme(plot.background=element_rect(fill="transparent"))


c <- ggplot(data,aes(x=Number,y=beta)) + geom_point()+ ggtitle(title)+ theme(plot.background=element_rect(fill="transparent"))
d <- ggplot(data,aes(x=Number,y=delta)) + geom_point()+ ggtitle(title) + theme(plot.background=element_rect(fill="transparent"))

return(list(a,b,c,d))
  
}


# function to convert dendrograms to phylos
dend2phylo <- function(dend) as.phylo(as.hclust.dendrogram(dend))


# Function to get tdend and other tree output from fineSTRUCTURE xml files
####### Input: data of type string of the filename like '*tree.xml''tree.xml'

readTrees <- function(treefile){
    treexml<-xmlTreeParse(treefile) ## read the tree as xml format
    ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
    ## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
    ## now is a good time to remove them via:
    #     ttree$node.label<-NULL
    ## Will will instead remove "perfect" node labels
    ttree$node.label[ttree$node.label=="1"] <-""
    ## And reduce the amount of significant digits printed:
    ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
    tdend<-myapetodend(ttree,factor=1,tol=1) # convert to dendrogram format
    
    ttreeclear <- ttree
    ttreeclear$node.label <-""  # Remove ALL labels
    tdendclear <- myapetodend(ttreeclear,factor=1,tol=1)

return(list("treexml"=treexml,"ttree"=ttree,"tdend"=tdend,"ttreeclear"=ttreeclear,"tdendclear"=tdendclear))
}



# Function to get boundaries of fineSTRUCTURE splits
####### Input: data of type 'new.split.matrix.A3'

getSplitBoundaries <- function(x,split,split_matrix){
    t <- split_matrix[x,split]-split_matrix[x+1,split]
    #print(t)
    if (t!=0) return(x) else return(NA)
} 


# Function to get XY coordinates of points on maps
####### Input: data of type 'split.id$map.name'. 
#Note it also needs the spatial polygon file 'espMapTotal' in the workspace to run properly.

getXYPosition <- function(loc,regloc) {
    print(loc)
    if ((is.na(loc))|(!loc%in%regloc$map.name)) {return(0)} else {
    test <- NA
    while (is.na(test[1])==T) {
    xpos <- regloc$V1[regloc$map.name==loc]+ runif(1,0,0.4)*regloc$area[regloc$map.name==loc]*sample(c(-1,1),1)
    ypos <- regloc$V2[regloc$map.name==loc]+ runif(1,0,0.4)*regloc$area[regloc$map.name==loc]*sample(c(-1,1),1)
  #test if in polygon overlaps with point
  test <- over(SpatialPoints(t(c(xpos,ypos)),proj4string=CRS(proj4string(espMapTotal))),espMapTotal[espMapTotal$map.name==loc,])
  }
    return(c(xpos,ypos))
  }
} 
 
## the same but for Spain B - more options (type) to allow for municiple, and provincial geo-codes
getXYPositionB <- function(loc,type,jitterSize=0) {
  #print(loc)
    if ((is.na(loc)==T) | (loc=='character(0)')) {return(0)} else {
    test <- NA
    while (is.na(test[1])==T) {
      
          if (type=="com") {
              xpos <- regloc$V1[regloc$map.name==loc]+ runif(1,0,0.4)*regloc$area[regloc$map.name==loc]*sample(c(-1,1),1)
              ypos <- regloc$V2[regloc$map.name==loc]+ runif(1,0,0.4)*regloc$area[regloc$map.name==loc]*sample(c(-1,1),1)
                    
              #test if in polygon overlaps with point
              test <- over(SpatialPoints(t(c(xpos,ypos)),proj4string=CRS(proj4string(espMapTotal))),espMapTotal[espMapTotal$map.name==loc,])
            }
          
          if (type=="prov") {
            
            #xpos <- jitter(provloc$V1[provloc$map.name==loc])
            #ypos <- jitter(provloc$V2[provloc$map.name==loc])
            xpos <- provloc$V1[provloc$map.name==loc]+ runif(1,0,0.4)*provloc$area[provloc$map.name==loc]*sample(c(-1,1),1)
            ypos <- provloc$V2[provloc$map.name==loc]+ runif(1,0,0.4)*provloc$area[provloc$map.name==loc]*sample(c(-1,1),1)
             
                    
              #test if in polygon overlaps with point
              test <- over(SpatialPoints(t(c(xpos,ypos)),proj4string=CRS(proj4string(espProvMap))),espProvMap[names(espProvMap)==loc])
            #test <- 1
            }
            if (type=="muni") {
                    if(!loc%in%muniloc$map.name){ return(0); test=1} else {
                        if(jitterSize>0) {
                          xpos <- jitter(muniloc$V1[muniloc$map.name==loc],amount=jitterSize)
                          ypos <- jitter(muniloc$V2[muniloc$map.name==loc],amount=jitterSize)
                            
                      #test if in polygon overlaps with point
                          test <- over(SpatialPoints(t(c(xpos,ypos)),proj4string=CRS(proj4string(espMuniMap))),espMuniMap) 
                          } else {
                              test <- 1
                              xpos <- muniloc$V1[muniloc$map.name==loc]
                              ypos <- muniloc$V2[muniloc$map.name==loc]
                          }   
                    }
              } 
      }
    return(c(xpos,ypos))
  }
} 

# For creating Spanish maps (not a function, but useful)
####### Input: espMapTotal, espComMap. 
 
espMapTotal$dummy <- 0
espProvTotal$dummy <- 0
espMuniTotal$dummy <- 0

  #myProjection = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  centreOfMadrid = espMapTotal["Madrid (Comunidad de)",]@polygons[[1]]@labpt
  myProjection = "+proj=laea +lat_0=40.5 +lon_0=-3.7 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  # project from centre of Madrid! 
  
  origProjection = proj4string(espMapTotal)
  espMapTotalProj = spTransform(espMapTotal,CRS(myProjection))
  espProvTotalProj = spTransform(espProvTotal,CRS(myProjection))
  espMuniTotalProj = spTransform(espMuniTotal,CRS(myProjection))
  ALLProj = spTransform(ALL,CRS(myProjection))

provLabels = as.character(rownames(espProvTotal@data))
provLabels[which(provLabels=="-lava")] <- "Álava"
provLabels[which(provLabels=="-vila")] <- "Ávila"
provLabelsProj = as.character(rownames(espProvTotalProj@data))
provLabelsProj[which(provLabelsProj=="-lava")] <- "Álava"
provLabelsProj[which(provLabelsProj=="-vila")] <- "Ávila"

canaryMelillaCom  <- which(rownames(espMapTotal@data)%in%c("Canarias","Melilla"))
canaryMelillaProv <- which(rownames(espProvTotal@data)%in%c("Santa Cruz de Tenerife","Palmas (Las)","Melilla"))

  commaplabels <- list("sp.text", coordinates(espMapTotal),as.character(espMapTotal@data$com.name,cex=0.6),col="black",alpha=0.5)
  commapborders <- spplot(espMapTotal,zcol="dummy",col.regions="transparent",colorkey=FALSE,main="title",col=hsv(0,0,0.2))
  commapborders2 <- spplot(espMapTotal[-canaryMelillaCom,],zcol="dummy",col.regions="transparent",colorkey=FALSE,main="title",col=hsv(0,0,0.2))
  
  commaplabelsProj <- list("sp.text", coordinates(espMapTotalProj),as.character(espMapTotalProj@data$com.name,cex=0.6),col="black",alpha=0.5)
  commapborders2Proj <- spplot(espMapTotalProj[-canaryMelillaCom,],zcol="dummy",col.regions="transparent",colorkey=FALSE,main="title",col=hsv(0,0,0.2))
  
provmaplabels <- list("sp.text", coordinates(espProvTotal),as.character(provLabels,cex=0.6),col="black",alpha=0.5)
  provmapborders <- spplot(espProvTotal,zcol="dummy",col.regions="transparent",colorkey=FALSE,main="title",alpha.regions=0.5,col=hsv(0,0,0.6))
  provmapborders2 <- spplot(espProvTotal[-canaryMelillaProv,],zcol="dummy",col.regions="transparent",colorkey=FALSE,main="title",alpha.regions=0.5,col=hsv(0,0,0.6))
  
  provmaplabelsProj <- list("sp.text", coordinates(espProvTotalProj),as.character(provLabels,cex=0.6),col="black",alpha=0.5)
  provmapborders2Proj <- spplot(espProvTotalProj[-canaryMelillaProv,],zcol="dummy",col.regions="transparent",colorkey=FALSE,main="title",alpha.regions=0.5,col=hsv(0,0,0.6))

 # munimapborders <- spplot(espMuniTotal2,zcol="dummy",col.regions="tranparent",alpha.regions=1,colorkey=FALSE,main="title",lwd=0.8)

a <- legend.master[order(legend.master$com.name),] #Because spplot assigns colors by order of the zcol, which is "com.name"
row.names(a) <- c()
 mapcolors1 <- a$colors[a$map.name%in%espMapTotal@data$map.name]
mapcolors2 <- rep("lightgrey",length(espMapTotal$com.name))
 
 canarias <-  sapply(espComMap@polygons, function(x)x@ID %in% c("Canarias"))
peninsulaPols <- espComMap[!canarias]
islandPols <- espComMap[canarias]
c <- 0.2
xmin=bbox(islandPols)[1,1]-c
ymin=bbox(islandPols)[2,1]-c
xmax=bbox(islandPols)[1,2]+c
ymax=bbox(islandPols)[2,2]+c

#canary_border <- layer(panel.rect(xleft=xmin,ybottom=ymin,
#                                    xright=xmax,
#                                    ytop=ymax))
 
canary_border <- list("sp.polygons",SpatialPolygons(list(Polygons(list(Polygon(cbind(c(xmin,xmin,xmax,xmax,xmin),
                                                                       c(ymin,ymax,ymax,ymin,ymin)))),"canary"))))

xylims <- list("Spain"=c(-13.1,5,34,44)) 
xylims3 <- list("Spain"=c(-9.5,5,36,44)) 
bboxOrig = makeBboxPolygon(espMapTotal,box = cbind(xylims3$Spain[c(1,3)],xylims3$Spain[c(2,4)]) );
bboxProj <- spTransform(bboxOrig,CRS(myProjection))
p=-10000
xylimsProj <- list("Spain"=c(bbox(bboxProj)[1,],bbox(bboxProj)[2,])+c(p,-30000,-20000,-p) )

xylims2 <- rbind(xylims$Spain[c(1,3)],xylims$Spain[c(2,3)],xylims$Spain[c(2,4)],
                 xylims$Spain[c(1,4)],xylims$Spain[c(1,3)])
bboxPolygon <- Polygons(list(Polygon(xylims2)),ID="boundbox")
bboxSpPoly <- SpatialPolygons(Srl=list(bboxPolygon),proj4string=CRS(proj4string(ALL)))
ALL_in_Spain <- which(over(ALL,bboxSpPoly)==1)

regioncol <- hsv(0,0,0.85)
Key <- as.data.frame("col")

#ALL_context_map <- spplot(ALL,zcol="dummy",xlim=subbox[1,],ylim=subbox[2,],colorkey=FALSE,col.regions= hsv(0,0,0.85),lwd=0.5,col="gray")
SpainMap_under <- spplot(ALL[ALL_in_Spain,],zcol="dummy",xlim=xylims$Spain[1:2],ylim=xylims$Spain[3:4],colorkey=FALSE,col.regions=regioncol,
                         par.settings = list(axis.line=list(col=NA)),sp.layout=(canary_border),main="title")
SpainMap_under2 <- spplot(ALL[ALL_in_Spain[-c(1,5)],],zcol="dummy",xlim=xylims3$Spain[c(1,2)],ylim=xylims3$Spain[c(3,4)],colorkey=FALSE,col.regions=regioncol,
                         par.settings = list(axis.line=list(col=NA)),main="title")

SpainMap_under2_sansSpain <- spplot(ALL[ALL_in_Spain[-7],],zcol="dummy",xlim=xylims3$Spain[c(1,2)],ylim=xylims3$Spain[c(3,4)],colorkey=FALSE,col.regions=regioncol,
                         par.settings = list(axis.line=list(col=NA)),main="title")


UnderMap <- SpainMap_under + layer(sp.polygons(islandPols,fill="lightgray"))
UnderMap2 <- SpainMap_under2
PortugalMap = ALL[ALL@data$NAME=="Portugal",]

# projected versions
f = ALLProj[ALL_in_Spain[-c(1,5)],]
s = ALLProj["222",]
# exclude Melilla from map
#mel = espMapTotalProj["Melilla",]; m = mel@polygons[[1]]@labpt 
#sapply(s@polygons[[1]]@Polygons,function(x) x@labpt)
f@polygons[[5]]@Polygons=s@polygons[[1]]@Polygons[-30]

SpainMap_under2Proj <- spplot(f,zcol="dummy",xlim=xylimsProj$Spain[c(1,2)],ylim=xylimsProj$Spain[c(3,4)],colorkey=FALSE,col.regions=regioncol,
                         par.settings = list(axis.line=list(col=NA)),main="title")

UnderMap2Proj <- SpainMap_under2Proj
PortugalMapProj = spTransform(PortugalMap,CRS(myProjection))

SpainPortugalMapProj = f
SpainWithoutPortMapProj = f[-4,]

SpainPortugalMapUnionProj = gUnaryUnion(SpainPortugalMapProj[-c(1,2),])
SpainMap_under2_portSpainUnion_Proj <- spplot(SpainPortugalMapUnionProj,xlim=xylimsProj$Spain[c(1,2)],ylim=xylimsProj$Spain[c(3,4)],colorkey=FALSE,col.regions=regioncol,
                         par.settings = list(axis.line=list(col=NA)),main="title")


#####################

# Function to get new split numbers for mapping
####### Input: data of type 'split.matrix.A3'. 

# Function to get an hclust object from fineSTRUCTURE tree (used in next function)

as.hclust.mine <- function (x, ...) 
{
   # if (!is.ultrametric(x)) 
      #  stop("the tree is not ultrametric")
    if (!is.binary.tree(x)) 
        stop("the tree is not binary")
    if (!is.rooted(x)) 
        stop("the tree is not rooted")
    n <- length(x$tip.label)
    x$node.label <- NULL
    bt <- sort(branching.times(x))
    inode <- as.numeric(names(bt))
    N <- n - 1L
    nm <- numeric(N + n)
    nm[inode] <- 1:N
    merge <- matrix(NA, N, 2)
    for (i in 1:N) {
        ind <- which(x$edge[, 1] == inode[i])
        for (k in 1:2) {
            tmp <- x$edge[ind[k], 2]
            merge[i, k] <- if (tmp <= n) 
                -tmp
            else nm[tmp]
        }
    }
    names(bt) <- NULL
    obj <- list(merge = merge, height = bt, order = 1:n, labels = x$tip.label, 
        call = match.call(), method = "unknown")
    class(obj) <- "hclust"
    obj
}

# creates matrix with splits in order of tree, and at each binary split the smallest split is given a new number
# maxsplits = integer
# tree = ttree object from the readTrees() function

getNewSplitNumbers <- function(maxsplits,tree){

if(class(tree)=="phylo") t <- as.hclust.mine(tree) else t <- tree

splits <- lapply(c(1:maxsplits),FUN=function(i) dendextend::cutree(t,i,use_labels_not_values=F,order_clusters_as_data=F,sort_cluster_numbers=TRUE))
 
 input <- as.data.frame(splits)

n <- nrow(input)
  j <- ncol(input)
  j1 <- j+1
  
  # Get split.points - list of indicies where each split occurs
 
    getSplitPoints <- function(x) {
          x1 <- x-1
          min(which(input[,x1]-input[,x]==-1))
          #rle(input[,x1])
    }

 split.points <- sapply(c(2:j),FUN=getSplitPoints)
 split.points <- c(0,split.points,n+1) # need to index either end of the list of samples (0, n+1)
 
 new.split.numbers <- input
 
  # generate new list of split codes
  
 for (t in c(2:j)) {
         tl <- t-1
   curr <-split.points[t]
   test.points <- split.points[c(1:tl,j1)]
   nearest <- test.points[which(abs(curr-test.points)==min(abs(curr-test.points)))]
   if (length(nearest) > 1) nearest <- nearest[1] # arbirarily choose the first one if two-equal
   direction <- sign(curr-nearest)

   new.split.numbers[,t] <- new.split.numbers[,tl]
    # Only change code of rows in smallest split group to value of t <==== Important!
         curr1 <- curr-1
         nearest1 <- nearest -1
    if (direction==1) new.split.numbers[c(nearest:curr1),t] <- t
    if (direction==-1) new.split.numbers[c(curr:nearest1),t] <- t
      
      print(c(t,curr,nearest,direction))
 
 }
names(new.split.numbers) <- paste("Split",seq(1,j),sep="_")

return(new.split.numbers)

}  

# Reduces the above tree split numbers so that splits < 5 individuals are ignored 
# This means that the total number of non-NA inds gets smaller as you go down the tree
getReducedSplitNumbers <- function(new.split.matrix){
  
  removeSmall <- function(x) {
        tab <- table(x)
        x2 <- x
        x2[which(as.character(x2)%in%names(tab)[which(tab < 5)])] <- NA
    return(x2)
  }
  new <- as.data.frame(apply(new.split.matrix[grep("Split",colnames(new.split.matrix))],MARGIN=2,FUN=removeSmall))
  new$ID <- new.split.matrix$ID
  return(new)
}


# Gets order of splits going down the tree (each binary pair)

getSplitOrder <- function(split.matrix) {
    splitOrder <- c(1,1)
    for (x in c(2:ncol(split.matrix))){
            splitI <- split.matrix[,x]
            splitI1 <- split.matrix[,(x-1)]
            
            this_split <- max(x)
            that_split <- unique(splitI1[which(splitI-splitI1>0)])
            splitOrder <- cbind(splitOrder,c(this_split,that_split))
     }
  return(splitOrder)
 }

## Exclude splits

excludeSplits <- function(split.matrix,exclusionSplits=NULL,excludeSamples=NULL){
if(!is.null(exclusionSplits)){
for (bad in exclusionSplits){
      exclusion_samples <- rownames(split.matrix[which(split.matrix[,bad]==bad),])
      split.matrix <- split.matrix[-which(rownames(split.matrix)%in%exclusion_samples),-bad]
      colnames(split.matrix) <- paste("Split_",c(1:ncol(split.matrix)),sep="")
    }
}

if(!is.null(excludeSamples)){
      split.matrix <- split.matrix[-which(rownames(split.matrix)%in%excludeSamples),]
      colsToRemove <- c()
      # if excluding these samples removes an entire split, remove the column too
      for(i in 1:ncol(split.matrix)){
        #print(i)
          if (length(which(split.matrix[,i]==i))==0) colsToRemove <- c(colsToRemove,i)
      }
             
      if (length(colsToRemove)>0) {
        print(paste("Removing splits: ",paste(colsToRemove,collapse=", ")))
        split.matrix <- split.matrix[,-colsToRemove]
      } else {print("No splits to remove :D")}
      
      colnames(split.matrix) <- paste("Split_",c(1:ncol(split.matrix)),sep="")
}
return(split.matrix)
}

load("~/Documents/ClareDPhil/DPhil/Spain/chromopainter/new.split.matrix.A2_v1_v7a.Rdata")
defaultBad=rownames(new.split.matrix.A2_v1_v7a)[new.split.matrix.A2_v1_v7a[,16]==16]


#######################
# Coloured dendrograms
#######################

SplitBoundaries<- function(j,split,split_matrix){
              #print(j)
              Vector <- split_matrix[,split]
              if(is.na(Vector[j])) return(NA) else {
                k=1
                t <- Vector[j]-Vector[j+k]
                while(is.na(t)){
                        t <- Vector[j]-Vector[j+k]
                        k=k+1
                }
                if (t!=0) return(j) else return(NA)
              }
}


getColouredDendro <- function(dend,m,split.matrix,fadeOutSplit=NULL,colourSizeFactor=15,collapse=F,colLevel=m,
                              rotate=F,removeClusters=NULL,zoom=NULL,flatten=F,returnDend=F,fadeOut=c(0.5,0.5),extra=0,extendLeaves=TRUE,...){
        .e <- environment()              
        
        if (!is.null(zoom)){
              stopifnot(class(dend)=="phylo")
              
              remove <- which(!rownames(split.matrix)%in%zoom)
              newTree =  drop.tip(dend,tip=rownames(split.matrix)[remove])
              dend = dendro_data(myapetodend(newTree,factor=1))
              
            }
        
      #  if (flatten){
          #    stopifnot(class(dend)=="phylo")
              
          #    dend = dendro_data(myapetodend(newTree,factor=1))
              
          #  }
                
        splitDend <- paste("Split",m,sep="_") # split level for the dendrogram
        
        split <- paste("Split",colLevel,sep="_") # split level for the colours
        
        # change segments info so that widths of leaves are equal
        if (collapse==T){
          
          if( class(dend)!="phylo") dend = dend2phylo(dend)
          
          if( length(table((split.matrix[,splitDend]))) > length(table((split.matrix[,split]))) ) splitHigher = splitDend else splitHigher=split
                       
          singles = sapply(unique(split.matrix[,splitHigher]),function(x) rownames(split.matrix)[which(split.matrix[,splitHigher]==x)][1])
          newTree =  drop.tip(dend,tip=dend$tip.label[!dend$tip.label%in%singles])
          dend = dendro_data(as.hclust.mine(newTree))
          
        }
        
        if(class(dend)=="phylo") {          
            r = myapetodend(dend,factor=1)
            treeOrder = labels(r)
        } else {
            treeOrder = as.character(dend$labels[,3])
        }
                  
        split.matrix <- split.matrix[match(treeOrder,rownames(split.matrix)),]          
        
        
        #split.matrix[is.na(split.matrix[,split]),split] <- 0
        split_order <- unique(split.matrix[,split])
            split_order <- split_order[!is.na(split_order)]
        boundaries <- sapply(1:(length(split.matrix[,split])-1),FUN=SplitBoundaries,split=split,split_matrix=split.matrix)
        split_positionsEnds <- c(boundaries[!is.na(boundaries)],nrow(split.matrix))
        split_positionsStarts <- c(1,split_positionsEnds[-length(split_positionsEnds)]+1)
        #split_positionsStarts <- c(0,split_positionsEnds[-length(split_positionsEnds)])
        split_cols <- numToColour[2,match(as.character(split_order),numToColour[1,])]
        if (!is.null(fadeOutSplit)) split_cols[which(!names(split_cols)%in%as.character(fadeOutSplit))] <- My.add.alpha(split_cols[which(!names(split_cols)%in%as.character(fadeOutSplit))],alpha=fadeOut[1],valueChange=fadeOut[2])
        split_positions <- as.data.frame(cbind(split_positionsStarts,split_positionsEnds,split_cols),stringsAsFactors=F)
        split_positions[,1] <- as.numeric(split_positions[,1])-0.5
        split_positions[,2] <- as.numeric(split_positions[,2])+0.5

        m = length(unique(split.matrix[,splitDend]))
        
        if(class(dend)=="phylo") dend = dendro_data(as.hclust.mine(dend))
                  
        if( (m*2-1) > length(which(dend$segments$y!=dend$segments$yend)) ) YMAX = 0 else YMAX <- sort(dend$segments$y[which(dend$segments$y!=dend$segments$yend)],decreasing=T)[(m*2-1)]
        
        # add root
        newdend <- dend$segments
        newdend <- rbind(newdend,c(mean(newdend$x[1],newdend$xend[1]),(newdend$y[1]+newdend$y[1]/50),mean(newdend$x[1],newdend$xend[1]),newdend$y[1]))
        newdend <- newdend[which(newdend$y>YMAX),]
                
        # collapse segments that are joined up at a higher part of the tree
        # not implemented        
        
        # extend lines to bottom of plot
        bottomBranches = newdend$yend<=YMAX
        
        if(extendLeaves) {
          
          newdend$yend[bottomBranches] <- 0
          
        } else {
          # cut everything off at the highest branch
          #newdend$yend[bottomBranches] <- max(newdend$yend[bottomBranches])
        }
        
        if(extendLeaves) newdend$yend[!bottomBranches] = newdend$yend[!bottomBranches] + extra  # lift only lower segments not part of bottom branches
        if(!extendLeaves) newdend$yend = newdend$yend + extra  # lift lower segments as well as bottom branches
        newdend$y = newdend$y + extra  # lift all upper segments up
        
        
        print(min(newdend$y))
        print(min(newdend$yend))
        
          
        newdendro = as.dendro(segments = newdend,labels = dend$labels,leaf_labels=dend$leaf_labels,class="dendro")
        
        
        
        # add root
      colour_bottom = -max(newdend$y)/colourSizeFactor
      nClusters <- nrow(split_positions)
        
      if(returnDend) {
        #newdend2 = as.dendro(segments = newdend,labels=)
        return(newdendro)
      } else {

      dendro_plot <- ggplot(newdend,environment=.e,rotate=T) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend),...) + theme_classic() +
                    geom_rect(data=split_positions,aes(ymin=colour_bottom,ymax=0,xmin=split_positionsStarts,xmax=split_positionsEnds),fill=split_positions$split_cols,colour=split_positions$split_cols) +
                    theme(axis.title=element_blank(),axis.text = element_blank(),
                    panel.background=element_blank(),
                    plot.background=element_blank(),
                    axis.line=element_blank(),
                    axis.ticks= element_blank())
       
        
        return(dendro_plot)
      }
}


######################
# remove blank pages from pdf
######################

removeBlanks <- function(filename,nPages){
  
  keepPages <- paste(seq(2,2*(nPages),2),collapse=" ")
  system(paste('mv ',filename,' temp2.pdf'))
  system(paste('pdftk temp2.pdf cat ',keepPages,' output ',filename,sep=""))
}


#######################
# Summary plots of chromo matrices with dendrograms
#######################

summaryMatrixPlotsAspect = 737/580
summaryMatrixPlots <- function(sum,split.matrix,tree,m,colourScale,title,plotDendro=T,legendCols=recieverLegend2) {
    
  if(plotDendro==T){
    g = getColouredDendro(dend=tree,m,split.matrix,collapse=T,rotate=T) +
              theme(axis.ticks.x=element_blank(),
              plot.margin = unit(c(0,0,0,0),"null"),
              panel.margin = unit(0,"null"),
              axis.ticks.length = unit(0,"null"),
              axis.ticks.margin = unit(0,"null")) +
              labs(x=NULL) + labs(y=NULL) + 
                        scale_x_continuous(expand=c(0,0)) +
                        scale_y_continuous(expand=c(0,0))
  }
     basep <- levelplot(raster(sum),col.regions=colourScale,colorkey=list(axis.text=list(cex=1.5)),cuts=50,scales=list(draw = FALSE),margin=F,
                              par.settings = theme.novpadding)
    
    # symbols for column labels
      shapes = rev(legendCols[rownames(sum),"shape"])
      cols = rev(legendCols[rownames(sum),"colour"])

    #### page parameters
         margin=0.03
         titleHeight = 0.09
         plotWidth=0.85
    
        myViewports <- list(viewport(0.5,0.5,height=(1-margin),width=(1-(margin/2))),
                            viewport(0.5,1,height=titleHeight,width=1,just="top"),
                            viewport(0.5,0,height=(1-titleHeight),just=c("bottom")),
                            viewport(1,0,height=1,width=plotWidth,just=c("right","bottom"))
                           )
    
    ##### start printing
         grid.newpage()
          pushViewport(myViewports[[1]])
          pushViewport(myViewports[[2]])
    
              grid.text(title,gp=gpar(cex=2))
         
          upViewport(1)
          pushViewport(myViewports[[3]])
          pushViewport(myViewports[[4]])

              print(basep,newpage=F)
              grid.points(x=rep(0.02,nrow(sum)),y=seq(1/(2*nrow(sum)),1,1/nrow(sum)),pch=shapes,size=unit(0.02,"npc"),gp=gpar(col=cols,fill=cols))
          
          upViewport(1)
            # grid.rect()
          pushViewport(viewport(x=0,0.5,
                        width=convertUnit(unit(1, "npc"), "npc",
                          "y", "dimension", "x", "dimension"),
                        height=convertUnit(unit((1-plotWidth), "npc"), "npc",
                          "x", "dimension", "y", "dimension"),
                        angle=90,just="top"))
            # grid.rect()
             print(g,newpage=F)
             
}

#######################
# Split-based box-plots
#######################

get_splits <- function(y,Split,split.matrix,indPairslist=indPairslist) {
      sp <- rownames(split.matrix[which(split.matrix[,Split]==y),])
      inds <- which((indPairslist$ID1%in%sp)&(indPairslist$ID2%in%sp))
      #print(y)
      return(inds)
}

get_pairs <- function(x,Split.matrix,indPairslist) {
    split <- paste("Split",x,sep="_")
    indPairslist$split <- 0
        
    l <- lapply(c(1:x),FUN=get_splits,Split=split,split.matrix=Split.matrix,indPairslist=indPairslist)
    
      for (j in c(1:x)) {  
      indPairslist[l[[j]],"split"] <- j
      }
    return(indPairslist$split)
}

SplitBasedBoxPlots <- function(filename,data,new.split.matrix,yMeasure,title=deparse(substitute(new.split.matrix)),type="pairs"){
    # Note: data must have ID rownames if type=NULL
  maxSplit <- max(as.numeric(gsub("[^[:digit:]]","",names(new.split.matrix))),na.rm=T)
   minSplit <- min(as.numeric(gsub("[^[:digit:]]","",names(new.split.matrix))),na.rm=T)
  if(type=="pairs"){
          indPairslist <- data[,c("ID1","ID2")] 
          pairssplitmatrix <- sapply(c(1:maxSplit),FUN=get_pairs,Split.matrix=new.split.matrix,indPairslist=indPairslist)
              colnames(pairssplitmatrix) <- paste("Split",seq(1,ncol(pairssplitmatrix)),sep="_")
          X <- cbind(data,pairssplitmatrix)
          X <- X[which(pairssplitmatrix[,maxSplit]!=0),]
    } else {X <- cbind(data,new.split.matrix[rownames(data),,drop=F]) }
    
    pdf(filename,width=12)
    for (splitNumber in c(minSplit:maxSplit)){
          split <- paste("Split",splitNumber,sep="_")
          X <- X[!is.na(X[,yMeasure]),]
          c <- X[,split]
          fillCols <- c(numToColour[2,match(as.character(unique(c)),numToColour[1,])],"transparent")
          names(fillCols) <- c(unique(c),"ALL")
          boxShapes <- c(as.numeric(unique(c))%%6,16)
          names(boxShapes) <- c(unique(c),"ALL")
          
          dummyvals <- c((splitNumber+1):maxSplit,"ALL")
          Xdummy <- X[1:length(dummyvals),]
          Xdummy[,yMeasure] <- NA
          Xdummy[,split]<- dummyvals

          X <- rbind(X,Xdummy)
          X[,split] <- factor(X[,split],levels=c(as.character(c(1:maxSplit)),"ALL"))
          X$endSplit <- "ALL"
          X$endSplit <- factor(X$endSplit, levels=levels(X[,split]))
          sample_sizes <- c(table(new.split.matrix[,split])[minSplit:splitNumber],rep(0,(maxSplit-splitNumber)),table(new.split.matrix[,"Split_1"]))
          names(sample_sizes) <- c(1:maxSplit,"ALL")
          
          g <- ggplot(X,aes_string(x=split,y=yMeasure,fill=split)) +
                    geom_boxplot()+
                    geom_boxplot(data=X,aes(x=endSplit,fill=endSplit))+
                      scale_fill_manual(values = fillCols,guide=F) +
                      #scale_shape_manual(values = boxShapes) +
                    scale_x_discrete(labels=sample_sizes,breaks=names(sample_sizes),name="split size") +
                    geom_hline(yintercept=median(X[,yMeasure],na.rm=T),linetype=3,col="red") +
                    geom_hline(yintercept=quantile(X[,yMeasure],na.rm=T,c(0.25,0.75)),linetype=3,col="black") +
                    ggtitle(title)
          print(g)
    }
dev.off()
}

#########################
# Functions to calculate tree differences

# first create cluster assigment strings for each split

splitsets<- function(split.matrix1,split.matrix2,s){
  n <- length(unique(split.matrix1[,s]))
  pw_diff <- matrix(nrow=n,ncol=n)
    
 for (i in c(1:n)){
   this_set1 <- rownames(split.matrix1)[which(split.matrix1[,s]==i)]
   for (j in c(1:n)){
        this_set2 <- rownames(split.matrix2)[which(split.matrix2[,s]==j)]
        overlap <- length(intersect(this_set1,this_set2))
        pw_diff[i,j] <- overlap
      }
  }
return(pw_diff)
}

#pw <- splitsets(split.matrix1=new.split.matrix.B_v2_v1,split.matrix2=new.split.matrix.B_v2_v2,10)


#########################
# Functions to plot splits on maps using pie-charts etc.

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

 ## Calculate contingency table for geography and splits, given a split matrix and a split level
get_contingency <- function(m,point.concordance,split.matrix,geog_level="map.name"){
 split <- paste("Split",m,sep="_")
 
  geogs <- point.concordance[match(rownames(split.matrix),point.concordance[,1]),geog_level]

  data <- as.data.frame(cbind(geogs,split.matrix[,split]),stringsAsFactors=FALSE)
  names(data) <- c("geography",split)
  data$geography[is.na(data$geography)] <- "na"
 
  tabulated <- as.matrix(table(data[,1],data[,2]))

  geogprops <- prop.table(tabulated,1)
  geogtotals <- margin.table(tabulated,1)
  splitprops <- prop.table(tabulated,2)
      splitprops_norm <- apply(tabulated,MARGIN=2,FUN=function(x) x/as.vector(geogtotals))
  splittotals <- margin.table(tabulated,2)
  
 return(list(m=m,tabulated=tabulated,geogprops=geogprops,geogtotals=geogtotals,splitprops=splitprops,
             splitprops_norm=splitprops_norm,splittotals=splittotals))
}
  
 ## Compute chi-squared terms for each cell (geog/split)
    
 #  chi <- (chisq.test(tabulated))
 # fish <- fisher.test(tabulated,simulate.p.value=TRUE)
  #  pearson <- chi$residuals
 #   pearson <- pearson[which(chi$expected<5)] # Ignore cells with expected count less than 5 (sampling too small)
 
# plotting function for proportions of each cluster in geographic region. Size of pie is sample size in that region
 
getPie <- function(x,rad,geog_data=regloc,geogtotals=geogtotals,
                            geogprops=geogprops,colors=colors,segOrder=NULL,coordNames=c("SPAIN.xcoords","SPAIN.ycoords")){
         geog_name <- x
      
      if(!geog_name%in%names(geogtotals)) {return(NA)} else {
       cx=0
       cy=0.02
       xy <- as.numeric(geog_data[which(geog_data$map.name==geog_name),coordNames])
          fig <- c()
          fig[1] <- xy[1]- rad + cx
          fig[2] <- fig[1] + 2*rad
          fig[3] <- xy[2]- rad + cy
          fig[4] <- fig[3] + 2*rad
         
        max <- max(geogtotals)
         total <- geogtotals[geog_name]
        radius = rad - (max-total)*0.06/max
       
       # Pie proportions
       topie <- geogprops[geog_name,]
       #colors[topie<0.01] <- "gray"
       if(is.null(segOrder)) segOrder <- c(1:length(topie))
       segOrder <- segOrder[!is.na(segOrder)]
       
       
       topie2 <- topie[segOrder]
       colors2 <- colors[match(segOrder,names(topie))]
       if(length(topie)==1) { topie2 <- topie; colors2 <- colors }
       
       return(list(fig=fig,radius=radius,colors2=colors2,topie2=topie2))
      }
  }


# returns and SpatialPolygons object
getPies2 <- function(y,centre,polyID,radius=0.5,adj=c(1,1),colours){
  
  edges = 200
  init.angle = 0
    
  x <- c(0, cumsum(y)/sum(y))
    dx <- diff(x)
    nx <- length(dx)
    twopi <- 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }
  polys <- list()
  for (i in 1L:nx) {
     # print(i)
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        a <- cbind(c(P$x, 0,P$x[1])*adj[1]+ centre[1,1], c(P$y, 0,P$y[1])*adj[2]+ centre[1,2])
        poly <- Polygons(list(Polygon(a)),ID=i)
        P <- t2xy(mean(x[i + 0:1]))
        #lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
        polys[[i]] <- poly
    }
  colours = colours[names(y)]
  
  spPoly <- SpatialPolygonsDataFrame(SpatialPolygons(polys,proj4string=CRS(proj4string(ALL))),data=data.frame(y),match.ID=F)
  plotOut <- spplot(spPoly,col.regions=colours,colorkey=F,par.settings = list(axis.line=list(col=NA)))
  
  spOut <- list("sp.polygons",SpatialPolygons(polys,proj4string=CRS(proj4string(ALL))),fill=colours,col="transparent",lwd=0)
  
  #spPoly <- list(Polygons(polys,polyID))
  # convert coordinate system
  #outPoly <- list("sp.polygons",spPoly)
  
  return(spOut)
}


plot.split.spread <- function(contingency_output=contingency,new.split.matrix,
                        point.concordance,pointsX="X",pointsY="Y",regionsmap=espMapTotal,baseMap,set="") {
      # Start plotting 
       # plot.new()
      if (set=="A"){
      can <- which(rownames(regionsmap@data)=="Canarias")
      regionsmap <- regionsmap[-can,]
        }
    
      splitprops_norm = contingency_output$splitprops_norm
      splitprops= contingency_output$splitprops
      splittotals= contingency_output$splittotals
      
      tabulated= contingency_output$tabulated
      m= contingency_output$m
      split <- paste("Split",m,sep="_")
      labels <- list("sp.text", coordinates(regionsmap),row.names(regionsmap@data),cex=0.6,col="black",alpha=0.5)
       
       plots <- list()
       treeorder <- match(row.names(new.split.matrix),point.concordance[,1])
        points <- SpatialPoints(point.concordance[treeorder,c(pointsX,pointsY)],proj4string=CRS(proj4string(regionsmap)))
        Mappoints <- list("sp.points",points,cex=0.8,lwd=1) 
              
       c <- as.character(new.split.matrix[,split])
       Mappoints$pch <- new.split.matrix[,split]%%6
       Mappoints$col <- numToColour[2,match(c,numToColour[1,])]
         mat <- match(row.names(regionsmap),rownames(splitprops_norm))
       
       for(split_include in c(1:m)) {
         if (!as.character(split_include)%in%colnames(splitprops_norm)) next # this split is too small to plot
         
        regionsmap$split_prop_norm <- 0
        regionsmap$split_prop_norm[!is.na(mat)] <- splitprops_norm[mat[!is.na(mat)],as.character(split_include)]
        regionsmap$split_prop_norm <- regionsmap$split_prop_norm/sum(regionsmap$split_prop_norm,na.rm=T) # Normalise to sum to one
        
       split.cols <- colorRampPalette(c("lightgray",My.add.alpha(numToColour[2,match(as.character(split_include),numToColour[1,])],valueChange=0.8,alpha=0.95)))(max(splittotals))
        mappoints <- Mappoints
       mappoints$col[c!=split_include] <- "darkgray"
        splitsize <- splittotals[as.character(split_include)]

        p<-spplot(regionsmap,alpha.regions=0.8,
                      zcol="split_prop_norm",col.regions=split.cols,at=c(seq(0,1,1/splittotals[as.character(split_include)]),1.000001), # normalise the scale
                      main=paste("cluster ",split_include," of ",m," (",splitsize,")",sep=""),colorkey=TRUE,xlim=xylims$Spain[1:2],ylim=xylims$Spain[3:4],lwd=0.5)
           
         P2 <- p + update(baseMap,col.regions="transparent") + update(p,sp.layout=list(labels,mappoints),col.regions="transparent",alpha=0)
         
         plot_name <- paste("plot",split_include,sep="_")
         plots[[plot_name]] <- P2
                 
      }
        return(plots)    
}
 
#contingency_output <- get_contingency(m=m,point.concordance.A,split.matrix)
plot.pies<- function(contingency_output=contingency,new.split.matrix=split.matrix,
                        point.concordance,pointsX="X",pointsY="Y",regionsmap=espMapTotal,
                        baseMap,pie_positions=regloc,multiple=T,set="",transparency=T) {
      # Start plotting 
       # plot.new()
      if (set=="A"){
      can <- which(rownames(regionsmap@data)=="Canarias")
      regionsmap <- regionsmap[-can,]
        }
    
      splitprops_norm = contingency_output$splitprops_norm
      splitprops= contingency_output$splitprops
      splittotals= contingency_output$splittotals
      
      tabulated= contingency_output$tabulated
      m= contingency_output$m
      split <- paste("Split",m,sep="_")
      labelPos <- coordinates(regionsmap)
      labelPos[,1] <- labelPos[,1]-0.5
      labels <- list("sp.text",labelPos,regionsmap@data$com.name,cex=0.8,col="transparent",alpha=0.5)

        regionsmap$dummy <- 0
        p <- spplot(regionsmap,zcol="dummy",colorkey=FALSE,col.regions="transparent",lwd=0.5,sp.layout=list(labels))
        forprint <- update(baseMap,aspect="fill",col.regions="lightgray",main=paste(m,"clusters")) + p

        colors <- numToColour[2,match(colnames(contingency_output$geogprops),numToColour[1,])]
        splitOrder <- getSplitOrder(new.split.matrix)
        toHighLight <- splitOrder[,m]
       # generate pie positions and sizes etc. for split m defined by contingency output
        pies <- lapply(pie_positions$map.name,FUN=getPie,geogprops=contingency_output$geogprops,
                       colors=colors,geogtotals=contingency_output$geogtotals,rad=0.07,geog_data=pie_positions,segOrder=as.character(unique(new.split.matrix[,m])))
        
       # get order of splits (for opacity)
      
       
       # Do one map with all colours
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(1,1)))
      print(forprint,vp = viewport(1,1))
      t2xy <- function(t,radius) {
        list(x = radius * cos(2*pi * t), y = radius * sin(2*pi * t))
      }
    
        for (this_pie in pies){ 
          if(is.na(this_pie[1])) next
           par(new=T,fig=this_pie$fig,mai=c(0,0,0,0),col=rgb(0,0,0,alpha=0.5),lty=0)
            #labels=c(pie_label,rep("",ncol(new.split.matrix)))
          this_pie$colors3 <- this_pie$colors2
          if((transparency==T)){
            
            #this_pie$colors3[which(!names(this_pie$topie2)%in%toHighLight)] <- add.alpha(this_pie$colors3[which(!names(this_pie$topie2)%in%toHighLight)],0.1)
                  this_pie$colors3[which(!names(this_pie$topie2)%in%toHighLight)] <- "darkgray"                       
                  P <- t2xy(seq.int(0, 1, length.out = 100),radius=this_pie$radius*15)
                  maxSeg <- max(this_pie$topie2/sum(this_pie$topie2))[1]
                  isMax <- this_pie$topie2[which(this_pie$topie2==1)]
                  if(length(isMax)==0) isMax <- 0
                  
                  if ((maxSeg==1)&(isMax%in%toHighLight)==T) {
                        pie(1,radius=this_pie$radius*15,edges=100,col=this_pie$colors3[which(this_pie$topie2==maxSeg)],xpd=NA,
                              cex=0.8,labels=NA,lty=0)                        
                        polygon(P$x,P$y,angle = 45,lty=1) 
                  } else {
                    pie(1,radius=this_pie$radius*15,edges=100,col="darkgray",xpd=NA,
                        cex=0.8,labels=NA,lty=0)
                    polygon(P$x,P$y,angle = 45,lty=1)
                    toColor <- which((names(this_pie$topie2)%in%toHighLight)&(this_pie$topie2!=0))
                    x <- this_pie$topie2
                    x <- c(0, cumsum(x)/sum(x))
                    dx <- diff(x)
                   
                        for (i in toColor) {
                        n <- max(2, floor(100 * dx[i]))
                        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),radius=this_pie$radius*15)
                        polygon(c(P$x, 0), c(P$y, 0),angle = 45,lty=1,col=this_pie$colors3[i])
                        }
                    
                  }

          }
      
      else {pie(this_pie$topie2,radius=this_pie$radius*15,edges=100,col=this_pie$colors3,xpd=NA,
              cex=0.8,labels=NA,lty=1)}
          
          pie_label <- as.character(regionsmap@data$com.name[match(names(this_pie$radius),rownames(regionsmap@data))])
          place <- t2xy(0.05,radius=this_pie$radius*15)
          text(1.3*place$x,1.3*place$y,pie_label,xpd=NA,adj=0)
          
       }
      
        if(multiple==T){
         # Start loop through split_includes
        for (split_include in c(1:m)){
           splitsize <- splittotals[as.character(split_include)]
           print(update(forprint,main=paste(split_include,"th cluster (",splitsize,")",sep="")))

        for (pie in pies){
          if(is.na(pie[1])) next
          colors3 <- pie$colors2
          colors3[which(!names(pie$topie2)%in%c(split_include))] <- "transparent"
      
       par(new=TRUE,fig=pie$fig,mai=c(0,0,0,0),cex=0.0001)
       pie(pie$topie2,radius=pie$radius*15,edges=50,labels="",col=colors3)
        }
     }
 }

}
 

# plotting function for showing tree splits as points on map

plotSplitsAsPoints <- function(filename=NULL,new.split.matrix,point.concordance,pointsX="X",pointsY="Y",
                               regionsmap=espMapTotal,printFile="yes",set="",baseMap,minSplitSize=5,dend=NULL,maxSplits=10){

  if (set=="A"){
      can <- which(rownames(regionsmap@data)=="Canarias")
      regionsmap <- regionsmap[-can,]
        }
  
  treeorder <- match(row.names(new.split.matrix),point.concordance[,1])
  points <- SpatialPoints(point.concordance[treeorder,c(pointsX,pointsY)],proj4string=CRS(proj4string(regionsmap)))
 nt <- min(length(grep("Split",names(new.split.matrix))),maxSplits)
 labels <- list("sp.text", coordinates(regionsmap),row.names(regionsmap@data),cex=0.6,col="black",alpha=0.5)
  
  plots <- list()
  if(printFile=="yes") pdf(filename,width=10) 
  
  regionsmap$dummy <- 0
  nPages=0
 for (m in c(1:nt)) {
   print(m)
    if(m>1){if(setequal(range(new.split.matrix[,m],na.rm=T),range(new.split.matrix[,m-1],na.rm=T))==TRUE) next} 
    # i.e if no change in non-NA splits then don't bother plotting.
    if(length(which(new.split.matrix[,m]==m))<minSplitSize) next #don't plot small clusters
    nPages=nPages + 1
    c <- as.character(new.split.matrix[,m])
     mappoints <- list("sp.points",points,col=numToColour[2,match(c,numToColour[1,])],pch=new.split.matrix[,m]%%6,cex=0.8) 
    
    this_map <- update(baseMap,main=paste(m,"clusters"),col.regions="lightgray") + 
        spplot(regionsmap,zcol="dummy",colorkey=FALSE,lwd=0.2,sp.layout=list(mappoints,labels),col.regions="transparent")
    
    if(!is.null(dend)) g <- getColouredDendro(dend,m,new.split.matrix)
   
    if(printFile=="yes"){
        grid.newpage()  
    print(this_map,vp=viewport(1,1))
    if(!is.null(dend)) print(g,vp = viewport(0.17,0.15,width=0.38,height=0.3),newpage=F)
        } else {
    plot_name <- paste("plot",m,sep="_")
    plots[[plot_name]] <- this_map
    
    }
 
   
   # print(g)
    }
  if(printFile=="yes"){
      dev.off()
      removeBlanks(filename,nPages)    
      
    } else {

  return(plots)

 }
}


makeScale <- function(values,colourSet,scalenum=100,nbreaks,fixedLims=NULL){
    par(mar=c(0,0,0,1.2))
    if(is.null(fixedLims)) colindex<-matrix(seq(min(values),max(values),length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
    if(!is.null(fixedLims)) colindex<-matrix(seq(fixedLims[1],fixedLims[2],length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
        
    colours <- colourScale(colindex[,1],colourSet)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=scalenum)
    scalephysicalpos <- seq(1,nbreaks,length.out=scalenum)
    image(1,1:nbreaks,colindex,xaxt="n",xlab="",yaxt="n",ylab="",col=colours,zlim=range(colindex),)
    axis(4,at=scalephysicalpos,labels=signif(scalelocs,2),las=2,cex.axis=0.5,tcl=-0.2,hadj=1)
    par(mar=opar$mar)
}

makeScaleWithHist <- function(values,colourSet,scalenum=100,nbreaks=70,fixedLims=NULL,binWidth=NULL,xlims=NULL,bgColour="transparent",xlab=NULL){
    par(mar=c(2,0,0,0),mgp=c(1,1,0))
    if(is.null(fixedLims)) colindex<-matrix(seq(min(values),max(values),length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
    if(!is.null(fixedLims)) colindex<-matrix(seq(fixedLims[1],fixedLims[2],length.out=scalenum),ncol=1,nrow=scalenum) # colour scale
   
    if(!is.null(binWidth)) getBreaks <- seq(min(values),max(values)+binWidth,by=binWidth) else getBreaks <- seq(min(values),max(values),length.out=nbreaks)
   
    h <- hist(values,breaks=getBreaks,plot=F)
    histMax <- max(h$density) + max(h$density)/20
    if(is.null(xlims)) xlims=c(min(values),max(values))
    
    hist(values,breaks=getBreaks,xlim=xlims,axes=F,xlab=xlab,ylab=NULL,main=NULL,probability=T,border="gray")
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
    axis(1,cex.axis=0.6,tck=-0.1,at=signif(seq(xlims[1],xlims[2],length.out=7),2))
     #box(lwd=0.5)    
    par(mar=opar$mar,mgp=opar$mgp)
}

## function to get adHocPoints for input into plotPontsOnMaps2 (i.e set morePoints=<output of getAdhocPoints>)

getAdhocPoints <- function(stuff,inds,m1,all=T,size=2){
  # if all==T provide all points in stuff$D but grayed out if not in inds
  # if all!=T provide only points in inds
  pointOrder <- c(which(!rownames(stuff$D)%in%inds),which(rownames(stuff$D)%in%inds))
  if (all==F) pointOrder <- which(rownames(stuff$D)%in%inds)
  inm1 <- which(rownames(stuff$D)[pointOrder]%in%inds)
  points <- SpatialPoints(stuff$D[pointOrder,c("Xmuni.mixture","Ymuni.mixture")],proj4string=CRS(proj4string(espMapTotal)))

  npoints <- nrow(points@coords)
  
  for(ms in m1){
   cols <- rep("darkgray",npoints)
      cols[inm1] <- recieverLegend2[m1,"colour"]
   shapes <- rep(16,npoints)
      shapes[inm1] <- recieverLegend2[m1,"shape"]
   bord <- rep("darkgray",npoints)
      bord[inm1] <- "black"
   sizes <- rep(1,npoints)
      sizes[inm1] <- size
  
  
   m1Points = list("sp.points",points,cex=sizes,col=bord,pch=shapes,fill=cols)
  return(m1Points)
}
}

## function to plot points on map of Spain (or any other map)

plotPointsOnMaps2 <- function(D,dend=NULL,zoom=F,split.matrix=NULL,m=1,
                              showLines=F,plotDendro=F,dendroPlot=NULL,baseMap=G,
                              levelsVector=NULL,levelSplits=10,
                              extraTitle="",pointSize=1,LWD=2,minSplitSize=0,
                              title=NULL,
                              trans=1,darken=1,showScale=F,colourSet=c("red","purple","blue"),
                              fixedLims=NULL,histXlims=NULL,vectorForHist=NULL,extraPoints=NULL,
                              mapLabels=T,newPage=T,numToColours=numToColour,numToShapes=NULL,pointLWD=0.2,
                              histBgCol="transparent",jitter="max",morePoints=NULL,scaleLab=NULL,
                              histBinWidth=0.03,bbox=NULL,pointsPCH=16,shapesVector=NULL,pointsBorderCol=NULL,defaultPointCol="black"){
                
                linesDummy <- sp.lines.tight[1:2,]
                linesDummy@lines[[1]]@Lines <- list(Line(cbind(c(-80,-80),c(-80,-80)))) # somewhere on the other side of the world
                linesDummy@lines[[2]]@Lines <- list(Line(cbind(c(-80,-80),c(-80,-80))))
                
                
                
                #Set points AND lines
                if (jitter=="max") {
                  pointColumns <- c("Xmuni.mixture","Ymuni.mixture") 
                  linesData <- sp.lines
                } 
                if (jitter=="tight") {
                  pointColumns <- c("Xmuni.ave.grand.tight","Ymuni.ave.grand.tight")
                  linesData <- sp.lines.tight
                }
                if (jitter==F) {
                  pointColumns <- c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")
                  linesData <- sp.lines.tight
                }
                if (jitter=="evenTighter") {
                  pointColumns <- c("Xmuni.ave.grand.evenTighter","Ymuni.ave.grand.evenTighter")
                  linesData <- sp.lines.tight
                }
                if (length(jitter)==2) {
                  pointColumns <- jitter 
                  linesData <- sp.lines
                } 
                
                if (zoom!=F) {
                  region = paste(gsub("\\.","; ",zoom),":",sep="")
                  dendroViewport = viewport(0.1,0.67,width=0.37,height=0.28,angle=90)
                  if (plotDendro==T) mapViewport = viewport(1,1,width=0.75,just=c("right","top"))
                  
                  # points <- SpatialPoints(D[,c("Xmuni.ave.grand.tight","Ymuni.ave.grand.tight")],proj4string=CRS(proj4string(espMapTotal)))
                   points <- SpatialPoints(D[,pointColumns],proj4string=CRS(proj4string(espMapTotal)))
                  if (showLines==T) linesToPlot <- linesData[which(rownames(linesData@data)%in%rownames(D)),] else linesToPlot=linesDummy
                  labels <- provmaplabels
                  #cods <- D[,c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")] 
                  #cods[which(is.na(cods[,1])),1] <- 0
                  #cods[which(is.na(cods[,2])),2] <- 0 
                  #points <- SpatialPoints(cods,proj4string=CRS(proj4string(espMapTotal)))
                  gal <- espMapTotal
                  if(zoom=="Galicia") gal@bbox <- espMapTotal["Galicia",]@bbox
                  if(zoom=="Basque")  gal@bbox <- espMapTotal[c("Rioja (La)","País Vasco","Navarra (Comunidad Foral de)"),]@bbox
                  if(zoom=="Cataluña.Aragón.Baleares")  gal@bbox <- espMapTotal[c("Aragón","Cataluña","Balears (Illes)"),]@bbox
                  if(zoom=="Andalucia.Murcia")  gal@bbox <- espMapTotal[c("Murcia (Región de)","Andalucía"),]@bbox
                  if(zoom=="other") gal@bbox <- bbox
                  
                  w <- gal@bbox[1,2]-gal@bbox[1,1]
                  h <- gal@bbox[2,2]-gal@bbox[2,1]
                  print(w/h)
                  if (plotDendro==F) mapViewport = viewport(0.5,0.5,width=1,height=w/h)
                 
                  baseMap <- update(baseMap,xlim=gal@bbox[1,]+c(-0.01,0.01),ylim=gal@bbox[2,]+c(-0.01,0.01))
                }
                if (zoom==F) {
                  region = ""
                  #dendroViewport = viewport(0.11,0.67,width=0.37,height=0.3,angle=90)
                  dendroViewport = viewport(1,0,width=0.37,height=0.3,just=c("right","bottom"))
                  mapViewport = viewport(0.5,1,just=c("centre","top"))
                  points <- SpatialPoints(D[,pointColumns],proj4string=CRS(proj4string(espMapTotal)))
                  
                  if (showLines==T) linesToPlot <- linesData[which(rownames(linesData@data)%in%rownames(D)),] else linesToPlot=linesDummy
                  labels=commaplabels
                }
                
                # is the projection of the points the same as baseMap?
                projBaseMap = proj4string(baseMap$panel.args.common$grid.polygons)
                if( projBaseMap != proj4string(points) ) {
                  print("warning: CRS of points doesn't match the baseMap. Updating points CRS.")
                  points = spTransform(points,CRSobj = CRS(projBaseMap) )
                  labels=commaplabelsProj
                }
                
                unknown <- D$Xmuni.ave.grand==0
                if (zoom==F) mappoints = list("sp.points",points,cex=rep(pointSize,nrow(D)),first=F)
                if (zoom!=F) mappoints = list("sp.points",points,cex=rep(pointSize*1.2,nrow(D)),first=F)
                #mappoints$cex[unknown] = 0.3
                mappoints$lwd = LWD
               
                ### print plots
   
                    if (!is.null(split.matrix)){
                          if ((length(grep("SPAIN",rownames(D)))>0)&(length(grep("SPAIN",rownames(split.matrix)))==0)) {
                            c <- split.matrix[match(rownames(D),paste("SPAIN_A2",rownames(split.matrix),sep="_")),as.numeric(m)]
                        } else {
                        c <- split.matrix[match(rownames(D),rownames(split.matrix)),as.numeric(m)]
                        }
                    } else {
                        c <- rep(1,nrow(D))
                    }
                        mappoints$col <- numToColours[2,match(as.character(c),numToColours[1,])]
                        print(c)
                        print(mappoints$col)
                         #mappoints$pch <- (c-1)%%6 + 15
                        if((is.null(numToShapes))&(is.null(shapesVector))) {
                          mappoints$pch <- numToShape[2,as.character(c%%dim(numToShape)[2])]
                        }
                        if(!is.null(shapesVector)) mappoints$pch <- shapesVector[match(rownames(D),names(shapesVector))] else mappoints$pch <- numToShapes[2,as.character(c)]                          
                        
                        print(table(mappoints$pch))
                        mappoints$fill <- mappoints$col
                        mappoints$lwd <- pointLWD
                        
                        #mappoints$bg <- mappoints$col
                        mappoints$col <- rep(defaultPointCol,length(c))
                        print(mappoints$col)
                
                        if(!is.null(pointsBorderCol)) mappoints$col <- pointsBorderCol[match(rownames(D),names(pointsBorderCol))]                          
                        print(mappoints$col)
                
                        tabc <- table(split.matrix[,m])
                        smallSplits <- names(tabc)[which(tabc<minSplitSize)]
                        mappoints$col[c%in%smallSplits] <- "transparent"
                        mappoints$fill[c%in%smallSplits] <- "transparent"
                        nSplitsShown <- length(tabc) - length(smallSplits)
                        
                       # print(mappoints)
                
                    if (!is.null(levelsVector)) {
                      print("making levelsVector colours")
                      if ( (length(grep("SPAIN_A2",names(levelsVector))) >0 ) & (length(grep("SPAIN_A2",rownames(D))) == 0) )  {
                        values <- levelsVector[match(paste("SPAIN_A2",rownames(D),sep="_"),names(levelsVector))]
                      } else {
                        values <- levelsVector[match(rownames(D),names(levelsVector))]
                      }
                        region = paste(extraTitle,region)
                        mappoints$col <- colourScale(values,colourSet=colourSet,nBreaks=100,fixedLims=fixedLims) # this must be the same as in the histogram
                        mappoints$pch <- pointsPCH
                        if( !is.null(shapesVector) ){
                          
                          mappoints$pch <- shapesVector[match(rownames(D),names(shapesVector))]
                          mappoints$fill <- mappoints$col
                          mappoints$col = defaultPointCol
                          
                        } else {
                          
                          if(pointsPCH > 20) {
                              mappoints$fill = mappoints$col 
                              mappoints$col = defaultPointCol
                              if(!is.null(pointsBorderCol)) mappoints$col <- pointsBorderCol[match(rownames(D),names(shapesVector))]
                          }
                      }
                        
                      }
    
                    if ((plotDendro==T)&(is.null(dendroPlot))) g <- getColouredDendro(dend,m,split.matrix,colourSizeFactor=20)
                    if (!is.null(dendroPlot)) g <- dendroPlot
                
                    mappoints$col[mappoints$col!="transparent"] <- My.add.alpha(mappoints$col[mappoints$col!="transparent"],alpha=trans,valueChange=darken)
                
                    if (is.null(title)) title = paste(region,m,"clusters")
                    if(is.null(morePoints)) morePoints <- list("sp.points",SpatialPoints(cbind(-80,-80),proj4string=CRS(proj4string(espMapTotal))))
                  
                
                    if (mapLabels==T) p <- update(baseMap, main = title, 
                                         colorkey = T) + spplot(linesToPlot, colorkey = T, col.regions = "darkgray", sp.layout = list(labels,morePoints, mappoints))
                    
                  #  print(title)
                  #  print(mapLabels)
                  #print(linetToPlot)
                    if (mapLabels==F) p <- update(baseMap, main = title, colorkey = T) + spplot(linesToPlot, colorkey = T, col.regions = "darkgray", sp.layout = list(morePoints,mappoints))
                    
                print(vectorForHist)
                if ( (!is.null(levelsVector))|(!is.null(vectorForHist)) ) record=F else record=T
                #graphics.off()
                plot.new()
                if (newPage==T) grid.newpage(recording=record)
                
                    if (showScale==T) {
                      mapViewport$height <- unit(0.90,"npc")
                      mapViewport$justification <- "top"
                      pushViewport(mapViewport)
                      print(p,newpage=F)
                      upViewport(1)
                      
                      #grid.rect()
                      #scaleViewport = viewport(1,0.5,width=0.06,height=0.9*p$aspect.ratio*mapViewport$height,just=c("right")
                      scaleViewport = viewport(0.05,0.01,height=0.1,width=unit(0.9*p$aspect.ratio,"npc"),just=c("left","bottom"))
                      pushViewport(scaleViewport)
                      par(new=T,fig=gridFIG())

                      
                      if (!is.null(levelsVector)&(!is.null(extraPoints))) values <- values[-match(extraPoints,names(values))]
                      if ( !is.null(vectorForHist)  ) values = vectorForHist
                        
                      print("making levelsVector histogram")
                      makeScaleWithHist(values,colourSet,scalenum=100,fixedLims=fixedLims,binWidth=histBinWidth,xlims=histXlims,bgColour=histBgCol,xlab=scaleLab)
                      
                    } else {
                      pushViewport(mapViewport)
                      print(p,newpage=F)
                    }
                    upViewport(1)
                    if (plotDendro==T) {
                      pushViewport(dendroViewport)
                      print(g,newpage=F)
                    }
                    
                    if (!is.null(extraPoints)){
      #                pushViewport(dendroViewport)
                       pushViewport(viewport(1,0,width=unit(1-0.9*p$aspect.ratio,"npc"),height=0.37,just=c("right","bottom")))
                       #grid.rect()
                       extras <- mappoints
                       extras$col[-which(rownames(D)%in%extraPoints)] <- "transparent"
                       extras$cex <- 0.5
                       p2 <- update(ALL_context_map,sp.layout=list(extras))
                       print(p2,newpage=F)
                    }
                    
                return()
        }



#######################
# Density plots (need an sp_grd)
#######################
espMapTotalSansCan <- espMapTotal[-which(rownames(espMapTotal@data)=="Canarias"),]
espMapTotalSansCanProj <- espMapTotalProj[-which(rownames(espMapTotalProj@data)=="Canarias"),]

plotDenisty <- function(density,colour="blue",title="you need a title",mainMap="SpainMap_under2",
                        mappoints=NULL,sp_grd,topAlpha=0.8,minDensity=0.00001,mapBounds=espMapTotalSansCan,notinthesea=NULL){
    # truncate density
    density[density<minDensity] <- 0
    sp_grd@data$dens <- density
    if(is.null(notinthesea)) notinthesea <- which(!is.na(over(sp_grd,mapBounds)[,1]))  # this is a bit slow, so only do once if have to
    hiDens <- which(sp_grd@data$dens>0)
    toPlot <- intersect(notinthesea,hiDens)
    if(length(toPlot)==0) {
      colour="transparent"
      print('no colour for this cluster')
      return(NULL)
    } else {
    
    if(length(colour)==1) colScale <- sapply(seq(0,topAlpha,length.out=500),FUN=function(x) My.add.alpha(colour,alpha=x))
    if(length(colour)>1) colScale=colorRampPalette(colour)(500)
    
    densityPlot <- spplot(sp_grd[toPlot,],zcol="dens",col.regions=colScale,colorkey=list(at=seq(0,1,length.out=50)),at=seq(0,1,length.out=500),col=NA,border=NA)
    
    if(!is.null(mappoints)) densityPlot <- update(densityPlot,sp.layout=list(mappoints))
    Map <- get(mainMap)
    
    pl <- update(Map,main=title) + densityPlot 
    
    if (mainMap=="SpainMap_under2") plot <- pl + provmapborders2 + commapborders2 else plot <- pl + provmapborders + commapborders
    
    rm(Map,pl,colScale,notinthesea)
    return(list(plot,densityPlot))
  }
}



plotDenistyContinuous <- function(density,colourSet=c("blue","purple","red"),
                        fixedLims=NULL,title="you need a title",mainMap="SpainMap_under2",
                        mappoints=NULL,sp_grd,topAlpha=0.8,minDensity=0.00001,
                        mapBounds=espMapTotalSansCan,notinthesea=NULL,transparencyVector=NULL){
    # truncate density
    density[density<minDensity] <- 0
    sp_grd@data$dens <- density
    if(is.null(notinthesea)) notinthesea <- which(!is.na(over(sp_grd,mapBounds)[,1]))
    hiDens <- which(sp_grd@data$dens>0)
    toPlot <- intersect(notinthesea,hiDens)
    
    if(is.null(fixedLims)) fixedLims = c(min(density),max(density))
    colorCuts =  seq(fixedLims[1],fixedLims[2],length.out=100)
    colorVector = colourScale(colorCuts,colourSet=colourSet,nBreaks=100,fixedLims=fixedLims)
    
    densityPlot = spplot(sp_grd[toPlot,],zcol="dens",col.regions=colorVector,cuts=100,at=colorCuts)
         
    Map <- get(mainMap)
    
    if(!is.null(transparencyVector)){      
      backgroundCol = Map$panel.args.common$col.regions
      # add alpha? ==> plot lightgray mask on top (or whatever the background colour on the map is)     
      sp_grd@data$trans <- transparencyVector
      transScale <- sapply(seq(topAlpha,0,length.out=500),FUN=function(x) My.add.alpha(backgroundCol,alpha=x))      
      
      trans = spplot(sp_grd[toPlot,],zcol="trans",col.regions=transScale,col=NA,border=NA)            
      densityPlot = densityPlot + trans
    }
                
    pl <- update(Map,main=title,col.regions="transparent") + densityPlot
    
    if(!is.null(mappoints)) pl <- update(pl,sp.layout=list(mappoints))
    
    if (mainMap=="SpainMap_under2") plot <- pl + provmapborders2 + commapborders2 else plot <- pl + provmapborders + commapborders
    
    rm(Map,notinthesea)
    return(list(plot,densityPlot))
}



###########
# a function to get split summary for density plots

getSplitSummary <- function(x,fun="sum",dens,split.matrix,splitColumn,weights=NULL){ 
  samples <- which(sampleIDs%in%rownames(split.matrix)[which(split.matrix[,splitColumn]==x)])
  if (length(samples)==0) summary <- rep(0,ncol(dens)) else {
  if (!is.null(weights)) dens <- dens*Weights
      if (fun=="sum")  {
        if (length(samples)==1) summary <- dens[samples,] else summary <- colSums(dens[samples,])
      }
      
      if (fun=="mean") {
        if (length(samples)==1) summary <- dens[samples,] else summary <- colMeans(dens[samples,])
      }
  }
  
    return(summary)
}

#######################
# Dendrogram stuff
#######################

setCols <- function(n,cols) {
        if(is.leaf(n)) {
          a <- attributes(n)
         # i <- i+1
          ind <- a$label
          attr(n, "nodePar") <- c(a$nodePar, list(lab.col = cols[ind]))
        }
        n
    }
#i <- 0

setColLab <- function(d,cols){
dendrapply(d,FUN=setCols,cols=cols)
}

### Remove leaves of a dendrogram  <---- Not working!

remLeaves <- function(n,IDs){
  if(is.leaf(n)) {
          a <- attributes(n)
          ind <- a$label
          if(!ind%in%IDs) {
           n
        }
    }
}

removeLeaves <- function(d,IDs){
    l <- dendrapply(d,FUN=remLeaves,IDs=IDs)
}

#colours <- numToColour[2,new.split.matrix.A_v5_v1$Split_5]
#names(colours)<- labels(TreeList_v5_v1$tdend)
#newdend <- setColLab(TreeList_v5_v1$tdend,cols=colours)

#plot(newdend)

#dend <- dendro_data(newdend)
#base <- ggplot(dend$segments) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + theme_bw()
#cur_split <- 5
#for (split in c(1:cur_split)){
 #   split_length <- new.split.matrix.A_v5_v1[,cur_split]
  #  base <- base + geom_segment(data=dend$segments[c(1:100),],aes(x=x, y=y, xend=xend, yend=yend),colour=) +

  
#plot.dendrogram(newdend)

#######################
# Sample QC stuff round(log(seq(0,0.2,0.001),10),2)
#######################

# A function to do standard analysis - input files of the form SPAIN.A.sample-stats

sampleQCplots <- function(inputfile,Lambda=10,title){
    
    
    S<- read.table(inputfile, header=T)
    S[,"logit_missing"]<- log(S[,"missing"])-log(1-S[,"missing"])
    
    Aberrant <- aberrant(S[,c("logit_missing","heterozygosity")],lambda=Lambda,uncorr=TRUE,standardize=FALSE)
    
    S$outliers <- 1
    S$outliers[Aberrant$outlier] <- 2
    h <- hist(S$heterozygosity,breaks=50,xlim=c(0.28,0.35),main=title)
      bigticks <- seq(-4, 0, by=1)
      smallticks <- log(c(seq(10^-4, 10^-3, by=10^-4),seq(10^-3, 10^-2, by=10^-3),seq(10^-2, 10^-1, by=10^-2)),10)
      labels <- sapply(bigticks, function(i) as.expression(bquote(10^ .(i))))
    m <- hist(log(S$missing,10),10,breaks=80,main=title,axes=F,xlab=("log10(missing)"),xlim=c(-4,0))
        axis(1,at=bigticks, labels=labels)
        axis(1,at=smallticks,labels=F)
        axis(2)
    p2 <- aberrant.plot(Aberrant)
    p1 <- plot(S$missing,S$heterozygosity,col=S$outliers,main=title)
          legend("bottomright",legend=c("All exclusions"),col="red",pch=1)
          abline(v=c(0.008,0.009,0.01,0.03))
    
    print(h)
    print(m)
    print(p1)
    print(p2)
}


#################################
# mixture modelling plots

donorPopnBars <- function(recieverPopn,mixmat,SUMMARYmat,barColour,title,popnLabels,...){
          par(las=1,...)
          x <- -mixmat[recieverPopn,]
          Labels <- sapply(names(x),getPopnLabel,popnLabels=popnLabels)[2,]
          max_raw <- max(SUMMARYmat/rowSums(SUMMARYmat))
          b <- barplot(x,xlim=c(-1,max_raw),col=barColour,horiz=T,axes=F,axisnames=F)
            segments(-percentile05bootstrap[recieverPopn,], b, -percentile95bootstrap[recieverPopn,],  b, lwd=2,xpd=NA)
            segments(-percentile05bootstrap[recieverPopn,], b - 0.1,  -percentile05bootstrap[recieverPopn,],b + 0.1, lwd=2,xpd=NA)
            segments(-percentile95bootstrap[recieverPopn,], b - 0.1,  -percentile95bootstrap[recieverPopn,],b + 0.1, lwd=2,xpd=NA)
            axis(1,at=c(-1,-0.5,0,0.5,max_raw),labels=c(1,0.5,0,0.5,round(max_raw,1)),outer=F)
            axis(2,at=b,labels=Labels,outer=F,pos=-0.95,lwd=0)
          par(new=T)
          barplot(SUMMARYmat[recieverPopn,]/sum(SUMMARYmat[recieverPopn,]),horiz=T,xlim=c(-1,max_raw),axes=F,axisnames=F,main=title)
}

receiverPopnBars <- function(donorPopn,receiversOnly=F,mixmat,SUMMARYmat,barColours,title=NULL,labels=F,text_labels=F,labelCex=1,...){
          innerMargin = 0.05
          margins <- c(0.7,2,0,0.5+innerMargin)
          par(las=1,mai=margins,...)
          if (receiversOnly==T) receivers <- grep("SpainPopn",rownames(mixmat)) else receivers=c(1:nrow(mixmat))
          x <- -mixmat[receivers,donorPopn]
          max_raw <- max(SUMMARYmat[receivers,]/rowSums(SUMMARYmat[receivers,]))
          Labels <- sapply(names(x),getPopnLabel,popnLabels=popnLabels)[2,]
          
          b <- barplot(x,xlim=c(-1,max_raw),col=barColours,horiz=T,axes=F,axisnames=F,xpd=NA,offset=-innerMargin)
            segments(-percentile05bootstrap[receivers,donorPopn]-innerMargin, b, -percentile95bootstrap[receivers,donorPopn]-innerMargin,  b,xpd=NA)
            segments(-percentile05bootstrap[receivers,donorPopn]-innerMargin, b - 0.1,  -percentile05bootstrap[receivers,donorPopn]-innerMargin,b + 0.1,xpd=NA)
            segments(-percentile95bootstrap[receivers,donorPopn]-innerMargin, b - 0.1,  -percentile95bootstrap[receivers,donorPopn]-innerMargin,b + 0.1,xpd=NA)
          
            #axis(1,at=c(-1,-0.5,0,0.5,1),labels=c(1,0.5,0,0.5,1),outer=F)
            axis(1,at=c(-1,-0.5,0)-innerMargin,labels=c(1,0.5,0),outer=F,xpd=NA)
          
            if (text_labels==T) axis(2,at=b,labels=Labels,outer=F,lwd=0)
            
          #par(new=T)
          #margins[2] <- innerMargin
          #par(mai=margins)
          par(new=T)
          barplot(SUMMARYmat[receivers,donorPopn]/rowSums(SUMMARYmat[receivers,]),horiz=T,xlim=c(-1,max_raw),axes=F,axisnames=F,main=title,offset=innerMargin)
          axis(1,at=c(0,0.5,max_raw)+innerMargin,labels=c(0,0.5,round(max_raw,1)),outer=F,xpd=NA)
          
          if (labels=="points") {
              split.numbers <- str_extract(rownames(mixmat)[receivers],"[[:digit:]]+")
              pointLabelsPch <- numToShape[2,as.character(as.numeric(split.numbers)%%dim(numToShape)[2])]
              pointLabelsCol <- numToColour[2,split.numbers]
              if (length(grep("donor",rownames(x)))>0){
                pointLabelsPch <- as.numeric(recieverLegend2[rownames(x),"shape"])
                pointLabelsCol <- recieverLegend2[rownames(x),"colour"]
              }
              points(y=b,x=rep(0,length(receivers)),pch=pointLabelsPch,col=pointLabelsCol,bg=pointLabelsCol,cex=labelCex,xpd=NA)
            }
 }


#################################
# function to filter samples for input into plotting functions (including filename string)

getFilteredPointsSpain <- function(maxDist=80000,version="close",zoom=F,split.matrix=NULL,colourVar="Clusters",filePrefix,excludePont=T,
                                   Portugal=F,addSPAINA2=T,projectPositions=FALSE){

      print(version)
      print(maxDist)
      print(zoom)
      
      if (is.null(split.matrix)) split.matrix <- point.concordance.A[,c(1:2)] # a dummy variable will all the samples
      
      if(version=="all") maxDist = NULL
      
        #if(version!="all") Distance_filter <- which(point.concordance.A[,paste(version,maxDist/1000,sep=".")]==0)
       if(version%in%c("close","far")) Distance_filter <- which(point.concordance.A[,paste(version,maxDist/1000,sep=".")]==0)
        
        if(version=="all") Distance_filter <- NULL
        if(version=="all.four.grand") Distance_filter <- which(is.na(point.concordance.A$maxDist))
        
    
         if(length(grep("SPAIN_A2",rownames(split.matrix)))>0) {
           QCexcluded <- which(!paste("SPAIN_A2_",rownames(point.concordance.A),sep="")%in%rownames(split.matrix))
         } else {
           QCexcluded <- which(!rownames(point.concordance.A)%in%rownames(split.matrix))
         }
         
         Pont <- c()
         if(excludePont==T) Pont <- which(((point.concordance.A$Xmuni.ave.grand.precise==point.concordance.A$Xmuni.ave.grand.precise[7])&(point.concordance.A$Ymuni.ave.grand.precise==point.concordance.A$Ymuni.ave.grand.precise[7])))
          
         toNotPlot <- unique(c(Distance_filter,QCexcluded,Pont))
         print(paste0("Excluding ",length(toNotPlot)," samples."))
          
         if (colourVar == "Clusters") filename = paste(filePrefix,'.splitsmap-grand.origin.',version,'-',maxDist/1000,'Km',sep="")
         if (colourVar == "Levels") filename = paste(filePrefix,'.levelsmap-',matrixName,'.grand.origin.',version,'-',maxDist/1000,'Km',sep="")
         
          if (zoom==F) {
            if(projectPositions) {
              print("Using projected versions")
              UnderMap2 = UnderMap2Proj
              commapborders2= commapborders2Proj
              provmapborders2 = provmapborders2Proj
            }
            baseMap = UnderMap2Proj + update(commapborders2Proj,col.regions=hsv(0,0,0.8),alpha.regions=0.5) + provmapborders2Proj + update(UnderMap2Proj,col.regions="transparent")
          }
          
          if (zoom!=F){
            gal <- espMapTotal
            if(zoom=="Galicia") gal@bbox <- espMapTotal["Galicia",]@bbox
            if(zoom=="Basque")  gal@bbox <- espMapTotal[c("Rioja (La)","País Vasco","Navarra (Comunidad Foral de)"),]@bbox
            if(zoom=="Cataluña.Aragón.Baleares")  gal@bbox <- espMapTotal[c("Aragón","Cataluña","Balears (Illes)"),]@bbox
            if(zoom=="Andalucia.Murcia")  gal@bbox <- espMapTotal[c("Murcia (Región de)","Andalucía"),]@bbox
            
            baseMap = spplot(gal,zcol="dummy",colorkey=F,main="title",par.settings = list(axis.line=list(col=NA))) + G
            filename = paste(filename,'.',zoom,'.pdf',sep="")
          }
        
      D <- point.concordance.A[-toNotPlot,]
      if(Portugal==T) {
       
        Port=point.concordance.PR_HM_NA[grep("POPRES",rownames(point.concordance.PR_HM_NA)),c("Person.ID","Xcountry_import","Ycountry_import")]
        colnames(Port) = c("Person.ID","Xmuni.mixture","Ymuni.mixture")
        D = rbind.fill(D,Port)
        rownames(D) = D$Person.ID
      }
      if(addSPAINA2){
        sp = !grepl("POPRES",rownames(D))
        rownames(D)[sp] = D$Person.ID[sp] = paste0("SPAIN_A2_",D$Person.ID[sp])
      }
      
      if(projectPositions){
          cols = grep("X|Y",colnames(D))
          newD = D
          for( i in cols[seq(1,length(cols),2)]){
            print(i)
              newPoints = spTransform(makeSpPoints(D[,c(i,i+1)],mapDat=espMapTotal),CRSobj = CRS(myProjection))
              newD[,c(i,i+1)] = newPoints@coords
          } 
          D = newD
      }
      
  return(list("D"=D,"filename"=filename,"baseMap"=baseMap))
         
}


#################################
# function to collapse small splits into their nearest 'parent' split, such that all splits are larger or equal to minSplitSize. 
# find nearest split using Tree - in phlyo format (note: this must be the tree corresponding to the split.matrix)

collapseSplits <- function(newSplits,minSplitSize,Tree){
    SplitsCollapsed <- c()
    tab <- table(newSplits)
    smallSplits <- names(tab)[which(tab<minSplitSize)]
    newSmallSplits <- smallSplits
    nSmall <- length(smallSplits)
    distances <- round(cophenetic.phylo(Tree),2)
    distMat <- matrix(NA,length(tab),length(tab))
        rownames(distMat) <- colnames(distMat) <- as.character(unique(newSplits))
    
    for (splita in rownames(distMat)){
      samplesa <- names(newSplits)[which(newSplits==as.numeric(splita))]
      
      for (splitb in rownames(distMat)[-c(1:which(rownames(distMat)==splita))]){
        samplesb <- names(newSplits)[which(newSplits==as.numeric(splitb))]
        dist <- distances[match(samplesa,rownames(distances)),match(samplesb,rownames(distances))]
      #  print(dist)
        if (length(samplesb==1)) value <- unique(dist) else value <- unique(unique(dist)[1,])
        
            if(length(value)!=1) {
              value <- value[1]
              #print(paste("Warning: more than one distance value between splits ",splita," and ",splitb,sep=""))
            }
        #print(paste(splita," and ",splitb,sep=""))
        distMat[match(splita,rownames(distMat)),match(splitb,rownames(distMat))] <- value
        distMat[match(splitb,rownames(distMat)),match(splita,rownames(distMat))] <- value
        }
    }
    
  while (length(smallSplits)>0) {
      for (split in smallSplits[order(names(tab[smallSplits]))]){
          SplitsCollapsed <- c(SplitsCollapsed,split)
          oldSmallSplits <- newSmallSplits
          #newSplit <- as.numeric(split)-1
          candidateSplits <- as.character(unique(newSplits))
          #smallerUniqueSplits <- uniqueSplits[which(uniqueSplits<=newSplit)]
         # print(candidateSplits)
          # find nearest split to merge with
          dist <- distMat[match(split,rownames(distMat)),match(candidateSplits,colnames(distMat))] # only look at splits that still exist!
          #print(dist)
          nearestSplit <- names(which(dist==min(dist,na.rm=T)))
          # merge with the smallest cluster
          candidates <- tab[nearestSplit]
          newSplitactual <- as.numeric(names(candidates)[which(candidates==min(candidates,na.rm=T))][1])
          if (length(candidates)>1) print(paste("more than one nearest cluster available for ",split,". Picking the smallest cluster: ",newSplitactual,sep=""))
          
          #newSplitactual=smallerUniqueSplits[which(abs(smallerUniqueSplits-newSplit)==min(abs(smallerUniqueSplits-newSplit)))]
 
          newSplits[which(newSplits==as.numeric(split))] <- newSplitactual
          
          print(paste("merging_",split,"_to_",newSplitactual,sep=""))
          tab <- table(newSplits)
          #print(tab)
          newSmallSplits <- names(tab)[which(tab<minSplitSize)]
          if (length(intersect(oldSmallSplits,newSmallSplits))==0) break  # don't continue collapsing splits as no need to
      }
      tab <- table(newSplits)
      smallSplits <- names(tab)[which(tab<minSplitSize)]
    }
  n <- length(SplitsCollapsed)
  print(paste("collapsed ",n," splits (out of ",nSmall," original small splits)",sep=""))
  print(SplitsCollapsed)
  return(newSplits)
}



collapseSplits2 <- function(newSplits,minSplitSize,Tree){
    SplitsCollapsed <- c()
    tab <- table(newSplits)
    smallSplits <- names(tab)[which(tab<minSplitSize)]
    newSmallSplits <- smallSplits
    nSmall <- length(smallSplits)
    
    
    
    distances <- round(cophenetic.phylo(Tree),2)
    distMat <- matrix(NA,length(tab),length(tab))
        rownames(distMat) <- colnames(distMat) <- as.character(unique(newSplits))
    
    for (splita in rownames(distMat)){
      samplesa <- names(newSplits)[which(newSplits==as.numeric(splita))]
      
      for (splitb in rownames(distMat)[-c(1:which(rownames(distMat)==splita))]){
        samplesb <- names(newSplits)[which(newSplits==as.numeric(splitb))]
        
        
        dist <- distances[match(samplesa,rownames(distances)),match(samplesb,rownames(distances))]
      #  print(dist)
        if (length(samplesb==1)) value <- unique(dist) else value <- unique(unique(dist)[1,])
        
            if(length(value)!=1) {
              value <- value[1]
              #print(paste("Warning: more than one distance value between splits ",splita," and ",splitb,sep=""))
            }
        #print(paste(splita," and ",splitb,sep=""))
        distMat[match(splita,rownames(distMat)),match(splitb,rownames(distMat))] <- value
        distMat[match(splitb,rownames(distMat)),match(splita,rownames(distMat))] <- value
        }
    }
    
  while (length(smallSplits)>0) {
      for (split in smallSplits[order(tab[smallSplits])]){
          SplitsCollapsed <- c(SplitsCollapsed,split)
          oldSmallSplits <- newSmallSplits
          #newSplit <- as.numeric(split)-1
          candidateSplits <- as.character(unique(newSplits))
          #smallerUniqueSplits <- uniqueSplits[which(uniqueSplits<=newSplit)]
         # print(candidateSplits)
          # find nearest split to merge with
          dist <- distMat[match(split,rownames(distMat)),match(candidateSplits,colnames(distMat))] # only look at splits that still exist!
          
          # get list of branch heights of MRCA between samples in split and all other samples
          canInds <- sapply(candidateSplits,function(x) names(newSplits)[newSplits==as.numeric(x)][1])
          
          mrcas <- sapply(canInds,function(x) {
            print(x)
            findMRCA(Tree, tips=c(names(newSplits)[newSplits==as.numeric(split)][1],x), type="node")
          })
          
          
          #print(dist)
          nearestSplit <- names(which(dist==min(dist,na.rm=T)))
          # merge with the smallest cluster
          candidates <- tab[nearestSplit]
          newSplitactual <- as.numeric(names(candidates)[which(candidates==min(candidates,na.rm=T))][1])
          if (length(candidates)>1) print(paste("more than one nearest cluster available for ",split,". Picking the smallest cluster: ",newSplitactual,sep=""))
          
          #newSplitactual=smallerUniqueSplits[which(abs(smallerUniqueSplits-newSplit)==min(abs(smallerUniqueSplits-newSplit)))]
 
          newSplits[which(newSplits==as.numeric(split))] <- newSplitactual
          
          print(paste("merging_",split,"_to_",newSplitactual,sep=""))
          tab <- table(newSplits)
          #print(tab)
          newSmallSplits <- names(tab)[which(tab<minSplitSize)]
          if (length(intersect(oldSmallSplits,newSmallSplits))==0) break  # don't continue collapsing splits as no need to
      }
      tab <- table(newSplits)
      smallSplits <- names(tab)[which(tab<minSplitSize)]
    }
  n <- length(SplitsCollapsed)
  print(paste("collapsed ",n," splits (out of ",nSmall," original small splits)",sep=""))
  print(SplitsCollapsed)
  return(newSplits)
}


#################################
# function to convert popn numbers to labels. Needs a list of popnLabels and pop should be of the for "donorPopn_x" or "SpainPopn_x"


getPopnLabel <- function(pop,popnLabels){
  
            names(popnLabels)[names(popnLabels)=="Neatherlands-Germany"] = "Netherlands-Sweden"
            names(popnLabels)[names(popnLabels)=="Spain_Galicia_coast"] = "Spain_Galicia_Pontevedra"
            names(popnLabels)[names(popnLabels)=="Spain_Galicia_inland"] = "Spain_Galicia_central"
            
            # Counts of individuals in donor group 27:
             # Belgium        Germany    Netherlands         Sweden  United Kingdom 
             #   1              6              9              8              1 
            
            Popn <- as.numeric(str_split(pop,"_",2)[[1]][2])
            type <- str_split(pop,"_",2)[[1]][1]
            pop2 <- NA
            if (type=="SpainPopn") {
                for (L in grep("Spain",names(popnLabels),value=T)){
                      if (Popn%in%popnLabels[[L]]) {
                        pop2 <- gsub("Spain_","",L)
                        break
                      }
                }
                
            }
            
              if (type=="donorPopn") {
                for (L in grep("Spain",names(popnLabels),value=T,invert=T)){
                      if (Popn%in%popnLabels[[L]]) {
                        pop2 <- L
                        break
                      }
                }

              }
    if(is.na(pop2)) {
      pop2 <- pop
      pop3 <- pop
    } else {
      pop3 <- paste(pop2,Popn,sep="_")
    }
        pop4 = pop3
            # some post-hoc adjustments for plotting!    
        pop2[grep("NorthAfrica.M-A-L",pop2)] = gsub("M-A-L","Mixed",pop2[grep("M-A-L",pop2)])
        pop3[grep("NorthAfrica.M-A-L",pop3)] = gsub("M-A-L","Mixed",pop3[grep("M-A-L",pop3)])   
        pop2[grep("Italy_5",pop3)] = paste0("Italy1")
        pop2[grep("Italy_12",pop3)] = paste0("Italy2")
        pop2[grep("Serbia_3",pop3)] = paste0("Serbia1")
        pop2[grep("Serbia_9",pop3)] = paste0("Serbia2")
        pop2[grep("Sub-saharan.Nigeria_80",pop3)] = paste0("Nigeria.YRI1")
        pop2[grep("Sub-saharan.Nigeria_43",pop3)] = paste0("Nigeria.YRI2")
        
        #pop2[grep("Sub-saharan.Nigeria",pop2)] = gsub("Sub-saharan.Nigeria","SS.Nigeria-YRI",pop2[grep("Sub-saharan.Nigeria",pop2)])
        #pop3[grep("Sub-saharan.Nigeria",pop3)] = gsub("Sub-saharan.Nigeria","SS.Nigeria-YRI",pop3[grep("Sub-saharan.Nigeria",pop3)])  
        
        pop3[grep("Sub-saharan.Nigeria",pop3)] = gsub("Sub-saharan.Nigeria","Nigeria.YRI",pop3[grep("Sub-saharan.Nigeria",pop3)])  
        
        pop2[grep("SS.Kenya-",pop2)] = gsub("SS.Kenya-","Kenya.",pop2[grep("SS.Kenya-",pop2)])
        pop3[grep("SS.Kenya-",pop3)] = gsub("SS.Kenya-","Kenya.",pop3[grep("SS.Kenya-",pop3)])  
        pop3[grep("SS.Kenya_76",pop3)] = gsub("SS.Kenya","Kenya",pop3[grep("SS.Kenya_76",pop3)])  
        
        
        #pop4[grep("Basque_318",pop4)] = "donorPopn_318"
        pop2[grep("Basque_318",pop4)] = "Basque1" # This is because the samples for 318 are labelled Basque1 in the Spain-only analysis :D
              
    return(c(pop2,pop3,pop4))
}

#list("NorthAfrica.M-A-L"="NAMIX","NorthMorocco"="NMOR","WesternSahara"="WESA")

getPopnLabel2 <- function(inds,split.matrix,popnLabels,type=1,prefix="SpainPopn_"){
    
    names(popnLabels)[names(popnLabels)=="Neatherlands-Germany"] = "Netherlands-Germany"
    names(popnLabels)[names(popnLabels)=="Spain_Galicia_coast"] = "Spain_Galicia_Pontevedra"
    names(popnLabels)[names(popnLabels)=="Spain_Galicia_inland"] = "Spain_Galicia_central"
            
    a= paste(prefix,split.matrix[inds,2],sep="")
    if(type==3) {
      out = a
    } else {
    out = sapply(a,getPopnLabel,popnLabels)[type,]
    }
    names(out) = inds    
    return(out)
  }


#################################
# function to Summaries of chunkcounts etc. given donor and recipient populations

getCountSummaries <- function(inputMatrix,DONOR.splits,SPAIN.splits,
         DONOR.split.matrix,SPAIN.split.matrix,
         DONOR.splitLevel,SPAIN.splitLevel,
         ALLpopns,ALLinds){
  
rowMax = rowSums(inputMatrix)
inputMatrix = inputMatrix[,rownames(DONOR.split.matrix)]

# re-normalise so all rows sum to the same (equivalent to normalising to sum to 1)
inputMatrix = rowMax*inputMatrix/rowSums(inputMatrix)

N <- length(ALLinds)

SUMMARYmat_1 <- matrix(nrow=N,ncol=length(DONOR.splits))
  colnames(SUMMARYmat_1) <- paste("donorPopn_",DONOR.splits,sep="")
  rownames(SUMMARYmat_1) <- ALLinds
SUMMARYmat_1mean <- SUMMARYmat_1

SUMMARYmat <- matrix(nrow=length(ALLpopns),ncol=length(DONOR.splits))
  colnames(SUMMARYmat) <- paste("donorPopn_",DONOR.splits,sep="")
  rownames(SUMMARYmat) <- ALLpopns
SUMMARYmat_mean <- SUMMARYmat

# populate summary matrix
print("getting info")
for (pop in DONOR.splits){
  inds <- rownames(DONOR.split.matrix)[which(DONOR.split.matrix[,DONOR.splitLevel]==pop)]
  inds <- inds[inds%in%colnames(inputMatrix)]
  mat <- inputMatrix[,inds]
  Pop <- paste("donorPopn_",pop,sep="")
  if (length(inds)==1) {SUMMARYmat_1[,Pop] <- mat} else {
    if (model=="sums"){
    SUMMARYmat_1[,Pop] <- rowSums(mat)  # NOTE: this is taking the sum (not mean) of chunk lengths (i.e total amount of genome copied from any donor from popn x)
    SUMMARYmat_1mean[,Pop] <- rowMeans(mat)
  }
  if (model=="means"){
    mat[mat==0] = NA
    SUMMARYmat_1[,Pop] <- rowMeans(mat,na.rm=TRUE)  # NOTE: this is taking the mean of chunk lengths (i.e average amount of genome copied from any donor from popn x)
    SUMMARYmat_1mean <- SUMMARYmat_1
  }
}
}

print('summarising')
for (pop in ALLpopns) {
  print(pop)
        Popn <- as.numeric(str_split(pop,"_",2)[[1]][2])
        type <- str_split(pop,"_",2)[[1]][1]
        if (type=="SpainPopn") {
          inds = rownames(SPAIN.split.matrix)[which(SPAIN.split.matrix[,SPAIN.splitLevel]==Popn)]
          if(length(grep("SPAIN_A2",inds))==0) inds[!grepl("POPRES",inds)] <- paste("SPAIN_A2",inds[!grepl("POPRES",inds)],sep="_")
      }
        if (type=="donorPopn") inds <- rownames(DONOR.split.matrix)[which(DONOR.split.matrix[,DONOR.splitLevel]==Popn)]
        inds <- inds[inds%in%rownames(SUMMARYmat_1)]
        mat <- SUMMARYmat_1[match(inds,rownames(SUMMARYmat_1)),]
        Meanmat <- SUMMARYmat_1mean[match(inds,rownames(SUMMARYmat_1mean)),]
        if(length(inds)>1) SUMMARYmat[pop,] <- colMeans(mat) else SUMMARYmat[pop,] <- mat
        if(length(inds)>1) SUMMARYmat_mean[pop,] <- colMeans(Meanmat) else SUMMARYmat_mean[pop,] <- Meanmat
}

return(list(SUMMARYmat_1,SUMMARYmat_1mean,SUMMARYmat,SUMMARYmat_mean))
}

# not finished!!!!
getCountSummaryDonor <- function(inputMatrix,DONOR.splits,DONOR.split.matrix){
  
donorSummaries <- list()
rowMax = rowSums(inputMatrix)
      
donorSummaries[[pop]]$predMat = predMat
donorSummaries[[pop]]$toFit = toFit

for(d in DONOR.splits){  
      pop = paste0("donorPopn_",d)      
      predInput = inputMatrix[,rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]!=d]] # exclude anyone from this donor population
      # renormalise
      predInput = rowMax*predInput/rowSums(predInput)
      
      predMat = sapply(DONOR.splits[DONOR.splits!=d],function(p){   
                  print(p)
                  inds <- rownames(DONOR.split.matrix)[which(DONOR.split.matrix[,DONOR.splitLevel]==p)]
                  inds <- inds[inds%in%colnames(inputMatrix)]
                  mat <- colMeans(sapply(DONOR.splits[DONOR.splits!=d],function(y){
                    inds2 = rownames(DONOR.split.matrix)[which(DONOR.split.matrix[,DONOR.splitLevel]==y)]
                    rowSums(predInput[inds,inds2])
                  }))
        })
        predMat = t(predMat)
        colnames(predMat) = rownames(predMat) = paste("donorPopn_",DONOR.splits[DONOR.splits!=d])
                  
      
  }
}

# version 2 for within-spain 
getCountSummaries2 <- function(inputMatrix,DONOR.splits,SPAIN.splits,
         DONOR.split.matrix,SPAIN.split.matrix,
         DONOR.splitLevel=2,SPAIN.splitLevel=2,
         ALLpopns,ALLinds,model="sums"){

  #NOTE: make sure the diagonals of inputMatrix are NA (if you have some individuals in rows and columns)
N <- length(ALLinds)

SUMMARYmat_1 <- matrix(nrow=N,ncol=length(DONOR.splits))
  colnames(SUMMARYmat_1) <- paste("SpainPopn_",DONOR.splits,sep="")
  rownames(SUMMARYmat_1) <- ALLinds
SUMMARYmat_1mean <- SUMMARYmat_1

SUMMARYmat <- matrix(nrow=length(ALLpopns),ncol=length(DONOR.splits))
  colnames(SUMMARYmat) <- paste("SpainPopn_",DONOR.splits,sep="")
  rownames(SUMMARYmat) <- ALLpopns
SUMMARYmat_mean <- SUMMARYmat

# populate summary matrix

for (pop in DONOR.splits){
  #print('computing summaries by individual')
  print(pop)
  inds1 <- rownames(DONOR.split.matrix)[which(DONOR.split.matrix[,DONOR.splitLevel]==pop)]
  inds = inds1[inds1%in%colnames(inputMatrix)]  
  if(length(inds1)>length(inds)) print('not all samples in DONOR matrix are in your input matrix')
  
  #print(head(inds))
  mat <- inputMatrix[,inds]
  Pop <- paste("SpainPopn_",pop,sep="")
  if (length(inds)==1) {SUMMARYmat_1[,Pop] <- mat} else {
    if (model=="sums"){
    SUMMARYmat_1[,Pop] <- rowSums(mat,na.rm=T)  # NOTE: this is taking the sum (not mean) of chunk lengths (i.e total amount of genome copied from any donor from popn x)
    SUMMARYmat_1mean[,Pop] <- rowMeans(mat,na.rm=T)
  }
  if (model=="means"){
    SUMMARYmat_1[,Pop] <- rowMeans(mat,na.rm=T)  # NOTE: this is taking the mean of chunk lengths (i.e average amount of genome copied from any donor from popn x)
    SUMMARYmat_1mean <- SUMMARYmat_1
  }
}
}

for (pop in ALLpopns) {
  #print('computing cross-popn summaries')
        Popn <- as.numeric(str_split(pop,"_",2)[[1]][2])
        type <- str_split(pop,"_",2)[[1]][1]
        inds <- rownames(SPAIN.split.matrix)[which(SPAIN.split.matrix[,SPAIN.splitLevel]==Popn)]
        if(sum(inds%in%colnames(inputMatrix))<length(inds)) print('not all samples in SPAIN matrix are in your input matrix')
        
        mat <- SUMMARYmat_1[rownames(SUMMARYmat_1)%in%inds,]
        Meanmat <- SUMMARYmat_1mean[rownames(SUMMARYmat_1mean)%in%inds,]
        if(length(inds)>1) SUMMARYmat[pop,] <- colMeans(mat) else SUMMARYmat[pop,] <- mat
        if(length(inds)>1) SUMMARYmat_mean[pop,] <- colMeans(Meanmat) else SUMMARYmat_mean[pop,] <- Meanmat
}

return(list(SUMMARYmat_1,SUMMARYmat_1mean,SUMMARYmat,SUMMARYmat_mean))
}

#################################
# function to plot data from mcmc uncertainty calculation
# Input:  output data from Garrett's code

plotClustUncertainty <- function(treeUnc,pointLocations,...){
  
  treeUncMatrix <- as.matrix(treeUnc[,-c(1,2)])
  uncertainSamples <- apply(treeUncMatrix,1,max)
  names(uncertainSamples) <-treeUnc[,1]
  plotPointsOnMaps2(pointLocations,levelsVector=uncertainSamples,showScale=T,newPage=T,...)
}


#################################
# function to rotate a matrix 90 degrees

matrot <- function(x) t(apply(x, 2, rev))


#################################
# function to plot heatmaps e.g for residuals
pngHeight = 10000

plotMixtureHeat <- function(toPlot,filename,title="Mixture proportions",cap=T,colourSet=NULL,Cex=3.5,ylab="recipients",
                            xlab="donors",negativeScale=T,nCols=100,fixedRange=NULL,plotPoints=T,recieverLegend2=recieverLegend2,
                            rowLabels=NULL,colLabels=NULL,pointScale=800,...){

toPlot2 <- toPlot

if((negativeScale==F)&(is.null(colourSet))) colourSet=c("white","yellow","red")
if(is.null(colourSet)) colourSet = c("red","white","blue")

if(cap==T){
  cap=mean(toPlot) + 2*sd(toPlot)
  print('capping matrix')
  toPlot2[toPlot2>cap] <- cap
  toPlot2[toPlot2<(-cap)] <- -cap
}

colourScale=colorRampPalette(colourSet,interpolate="linear")(nCols)
#cap=5

Range <- max(abs(toPlot2),na.rm=T)
if(negativeScale==F) Range=c(min(toPlot2,na.rm=T),max(toPlot2,na.rm=T))
if(negativeScale==T) Range=c(-Range,Range)
if(!is.null(fixedRange)) Range=fixedRange

if(is.null(rowLabels)) rowLabels <- sapply(rownames(toPlot2),getPopnLabel,popnLabels=popnLabels)[1,]
if(is.null(colLabels)){
  colLabels <- sapply(gsub(".null","",colnames(toPlot2)),getPopnLabel,popnLabels=popnLabels)[1,]
  colLabels[grep("null",colnames(toPlot2))] <- paste(colLabels[grep("null",colnames(toPlot2))],".null",sep="")
}

if(plotPoints==T){      
  rowShapes <- as.numeric(recieverLegend2[rownames(toPlot2),"shape"])
  rowColours <- recieverLegend2[rownames(toPlot2),"colour"]
 
  colShapes <- as.numeric(recieverLegend2[gsub(".null","",colnames(toPlot2)),"shape"])
  colColours <- recieverLegend2[gsub(".null","",colnames(toPlot2)),"colour"]
  }
png(filename,...)
      layout(mat=matrix(c(1,2),1,2),widths=c(0.9,0.1))
      par(mar=c(11,11,2,0.5),mgp=c(10,1,0),cex=Cex)
      image(as.matrix(toPlot2),col=colourScale,zlim=Range,axes=F,main=title,
            xlab=xlab,ylab=ylab)
      axis(1,at=seq(0,1,length.out=nrow(toPlot2)),labels=rowLabels,las=2,tick=F)
      axis(2,at=seq(0,1,length.out=ncol(toPlot2)),labels=colLabels,las=2,tick=F)
      if(plotPoints==T){
      points(x=seq(0,1,length.out=nrow(toPlot2)),y=rep(-nrow(toPlot2)/pointScale,nrow(toPlot2)),col="black",pch=rowShapes,xpd=NA,bg=rowColours)
      points(y=seq(0,1,length.out=ncol(toPlot2)),x=rep(-ncol(toPlot2)/pointScale,ncol(toPlot2)),col="black",pch=colShapes,xpd=NA,bg=colColours)
      }
      #1300,1933
      par(mar=c(11,0,2,4))
      colindex<-t(matrix(seq(Range[1],Range[2],length.out=nCols),ncol=1,nrow=nCols))
      image(1,1:nCols,colindex,xaxt="n",yaxt="n",xlab="",ylab="",col=colourScale,zlim=range(colindex))
      scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=10)
      axis(2,at=seq(1,nCols,length.out=10),labels=signif(scalelocs,2),las=2,cex.axis=1,side=4)

dev.off()
par(opar)
}


#################################
# calculate 2-norm of vector x
norm_vec <- function(x) sqrt(sum(x^2))

#################################
# Read in fineSTRUCTURE/CHROMOPAINTER data

readChromo = function(prefix,chromv,extra=""){
        ch <- paste("./chromo_out/",prefix,"_v",chromv,extra,".chunkcounts.out",sep="")
        l <-  paste("./chromo_out/",prefix,"_v",chromv,extra,".chunklengths.out",sep="")
        mut <- paste("./chromo_out/",prefix,"_v",chromv,extra,".mutationprobs.out",sep="")
          
        d <-fixNames(as.matrix(read.table(ch,row.names=1,header=T,skip=1))) # read in chunk counts etc.
        c<- fixNames(as.matrix(read.table(l,row.names=1,header=T))) 
        mt<- fixNames(as.matrix(read.table(mut,row.names=1,header=T))) 
        
        toFix = c(grep("150X$",rownames(d)),grep("150X$",colnames(d)))
        
            rownames(d) <- gsub("150X$","150XX",rownames(d))
            colnames(d) <- gsub("150X$","150XX",colnames(d))
            rownames(c) <- gsub("150X$","150XX",rownames(c))
            colnames(c) <- gsub("150X$","150XX",colnames(c))
            rownames(mt) <- gsub("150X$","150XX",rownames(mt))
            colnames(mt) <- gsub("150X$","150XX",colnames(mt))
          
        meanc <- c/d
        diag(meanc) <- 0
          
        assign(paste('chunkcountsmatrix_v',chromv,extra,sep=""),d)
        assign(paste('chunklengthsmatrix_v',chromv,extra,sep=""),c)
        assign(paste('meanchunklengthsmatrix_v',chromv,extra,sep=""),meanc)
        assign(paste('mutationprobsmatrix_v',chromv,extra,sep=""),mt)
          
        save(list=c(paste('chunkcountsmatrix_v',chromv,extra,sep=""),paste('chunklengthsmatrix_v',chromv,extra,sep=""),
                      paste('meanchunklengthsmatrix_v',chromv,extra,sep=""),
                      paste('mutationprobsmatrix_v',chromv,extra,sep="")),
               file=paste(prefix,'.ChromopainterOutput_v',chromv,extra,'.Rdata',sep=""))
        
        print(paste('Objects saved to ',prefix,'.ChromopainterOutput_v',chromv,extra,'.Rdata',sep=""))
}

readFine = function(prefix,chromv,finev,new=FALSE){
  m <- paste("./chromo_out/",prefix,"_v",chromv,"_v",finev,".mcmc.xml",sep="") ## finestructure mcmc file
  if(new) m <- paste("./chromo_out/",prefix,"_v",chromv,"_v",finev,".mcmc.new.xml",sep="") ## finestructure mcmc file  
    Mc1 <- xmlTreeParse(m)
    Mc2 <- as.data.frame.myres(Mc1)
    assign(paste('mcmcxml_v',chromv,'_v',finev,sep=""),Mc1)
    assign(paste('mcmcdata_v',chromv,'_v',finev,sep=""),Mc2)

print(paste('reading trees...'))

  t <- paste("./chromo_out/",prefix,"_v",chromv,"_v",finev,".tree.xml",sep="") ## finestructure tree file
  if(new) t <- paste("./chromo_out/",prefix,"_v",chromv,"_v",finev,".tree.new.xml",sep="") ## finestructure tree file
  
    Tr <- readTrees(t)
    assign(paste('TreeList_v',chromv,'_v',finev,sep=""),Tr)
    
    meancoi <- paste("./chromo_out/",prefix,"_v",chromv,"_v",finev,".meancoincidence.xml",sep="") # pairwise coincidence, .i.e. proportion of MCMC files where individuals are found in the same 
        pairwise <- as.matrix(read.csv(meancoi,row.names=1)) # read in the pairwise coincidence file we created earlier
    #copy <- paste("./chromo_out/",prefix,suffix2,".copying.xml",sep="") # population copying matrix
    #    copying <- as.matrix(read.csv(copy,row.names=1))
    #copynorm <- paste("./chromo_out/",prefix,suffix2,".copyingnorm.xml",sep="") # population copying matrix
    #    copyingnorm <- as.matrix(read.csv(copynorm,row.names=1))
        pairwise <- fixNames(pairwise)
  
    map <-extractValue(Tr$treexml,"Pop") # map state as a finestructure clustering
    maplist<-popAsList(map) # .. and as a list of individuals in populations
    fullorder<-labels(Tr$tdend) # the order according to the tree
    mapstatematrix<-groupingAsMatrix(maplist)[fullorder,fullorder] 

    assign(paste('mapstate_v',chromv,'_v',finev,sep=""),map)
    assign(paste('mapstatelist_v',chromv,'_v',finev,sep=""),maplist)    
    assign(paste('mapstatematrix_v',chromv,'_v',finev,sep=""),mapstatematrix)
    
    assign(paste('PairwiseCoincidence_v',chromv,'_v',finev,sep=""),pairwise)
  #  assign(paste('copying_v',chromv,'_v',finev,sep=""),copying)
   # assign(paste('copyingnorm_v',chromv,'_v',finev,sep=""),copyingnorm)
    
    if(new){
      save(list=c(paste('TreeList_v',chromv,'_v',finev,sep="")),file=paste(prefix,'.TreeList_v',chromv,'_v',finev,'.new.Rdata',sep=""))
      save(list=c(paste('mapstatelist_v',chromv,'_v',finev,sep=""),paste('mcmcdata_v',chromv,'_v',finev,sep=""),paste('mcmcxml_v',chromv,'_v',finev,sep=""),paste('PairwiseCoincidence_v',chromv,'_v',finev,sep="")),file=paste(prefix,'.McmcData_v',chromv,'_v',finev,'.new.Rdata',sep=""))
  
    }else {
        save(list=c(paste('TreeList_v',chromv,'_v',finev,sep="")),file=paste(prefix,'.TreeList_v',chromv,'_v',finev,'.Rdata',sep=""))
        save(list=c(paste('mapstatelist_v',chromv,'_v',finev,sep=""),paste('mcmcdata_v',chromv,'_v',finev,sep=""),paste('mcmcxml_v',chromv,'_v',finev,sep=""),paste('PairwiseCoincidence_v',chromv,'_v',finev,sep="")),file=paste(prefix,'.McmcData_v',chromv,'_v',finev,'.Rdata',sep=""))
    }
   
}
  
newSplitMatrix = function(pre,chromv,finev,prefix,maxsplitboundaries=20,new=FALSE){
  
if(new) load(paste(prefix,'.McmcData_v',chromv,'_v',finev,'.new.Rdata',sep="")) else  load(paste(prefix,'.McmcData_v',chromv,'_v',finev,'.Rdata',sep=""))

if(new) load(paste(prefix,'.TreeList_v',chromv,'_v',finev,'.new.Rdata',sep="")) else load(paste(prefix,'.TreeList_v',chromv,'_v',finev,'.Rdata',sep=""))

m <- paste('mapstatelist_v',chromv,'_v',finev,sep="")
  print(paste("finev",finev))
  mapstatelist <- get(m)
  TreeList <- get(paste('TreeList_v',chromv,'_v',finev,sep=""))
 nt <- length(mapstatelist)
 new.split.matrix <- getNewSplitNumbers(nt,TreeList$ttree)
  
     new.split.matrix.b <- getReducedSplitNumbers(new.split.matrix)
     nsplits <- min(ncol(new.split.matrix),maxsplitboundaries)
     s <- sapply(c(1:(nrow(new.split.matrix)-1)),FUN=getSplitBoundaries,split=paste("Split_",nsplits,sep=""),split_matrix=new.split.matrix)
     split_boundaries <- s[!is.na(s)]
      
  if(new){
       assign(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,sep=""),new.split.matrix)
       assign(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,'b',sep=""),new.split.matrix.b)
       assign(paste('split_boundaries.',pre,'_v',chromv,'_v',finev,sep=""),s)
       save(list=c(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,sep=""),paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,'b',sep="")),
              file=paste('./new.split.matrix.',pre,'_v',chromv,'_v',finev,'.new.Rdata',sep=""))
       save(list=paste('split_boundaries.',pre,'_v',chromv,'_v',finev,sep=""),
              file=paste('split_boundaries.',pre,'_v',chromv,'_v',finev,'.new.Rdata',sep=""))
} else {
     assign(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,sep=""),new.split.matrix)
     assign(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,'b',sep=""),new.split.matrix.b)
     assign(paste('split_boundaries.',pre,'_v',chromv,'_v',finev,sep=""),s)
     save(list=c(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,sep=""),paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,'b',sep="")),
            file=paste('./new.split.matrix.',pre,'_v',chromv,'_v',finev,'.Rdata',sep=""))
     save(list=paste('split_boundaries.',pre,'_v',chromv,'_v',finev,sep=""),
            file=paste('split_boundaries.',pre,'_v',chromv,'_v',finev,'.Rdata',sep=""))
  }
}

getlabels <- function(y,input) { 
  w <- which(input==y)[floor(length(which(input==y))/2)] 
  if (length(w)==0) {w <- which(input==y)}
  return(w)
}


#################################
# plot fineSTRUCTURE dendrograms and coancestry matrices

colouredDendrograms <- function(plotdir,prefix,pre,chromv,finev,some.colorsEnd2=NULL){
  
list <- get(paste('TreeList_v',chromv,'_v',finev,sep=""))
dend <- dendro_data(list$tdend)
split.matrix <- get(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,sep=""))
split.matrixb <- get(paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,'b',sep="")) 

pdf(paste(plotdir,'/',prefix,'.dendrograms.v',chromv,'_v',finev,'.pdf',sep=""),width=12)
  
split_cols <- numToColour[2,match(split.matrixb[,ncol(split.matrixb)],numToColour[1,])]
nsplits <- length(colnames(split.matrixb))-1
          for (i in c(1:nsplits)) {
             
            split <- paste("Split",i,sep="_")
            text <- as.character(split.matrix[,split])
            text[is.na(text)] <- "na"
           t <-  sapply(text,FUN=getlabels,input=text)
            text[-t[!is.na(t)]] <- ""
            
            split_cols <- numToColour[2,match(split.matrixb[,split],numToColour[1,])]
          
           g<- ggplot(dend$segments) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + theme_bw() +
                    scale_x_continuous(breaks=c(1:length(dend$labels$label)),labels=text) + 
                    theme(axis.text.y = element_text(color="white"),
                    axis.ticks.x = element_line(color=split_cols,size=2),
                    axis.ticks.y = element_line(color="white"),
                    axis.ticks.length=unit(20,"points"),axis.title=element_text(color="white"),
                    panel.background=element_rect(fill="white")) 
                    ggtitle(paste(split))
            print(g)
          }
      dev.off()
}



plotFineMatrix <- function(matrixTypes,chromv,plotdir,prefix,pre=NULL,finev=NULL,data=NULL,
                           tree=T,otherOrder=NULL,treeCut=7000,
                           plotOrderRows=NULL,plotOrderCols=NULL,
                           labelPoints=NULL,title="",labels=NULL,split_boundaries=NULL,cap=T,
                           some.colorsEnd2=NULL,picSize=c(2200,2000)){
  # get colours
    if(is.null(some.colorsEnd2)) some.colorsEnd2 <- colorRampPalette(c("yellow","red","blue","black"),interpolate="linear")(100) # as above, but with a dark grey final for capped values
  
  # get row/col order of raw matrices  
    if(is.null(data)){
              data1 <- get(paste(matrixTypes[1],'_v',chromv,sep=""))              
    } else {
      data1 =get(data)
    }
              rawdataRowOrder = rownames(data1)
              rawdataColOrder = colnames(data1)

  # get order for rows (columns in plot)
      if(tree==T){
        print('getting tree info')
          treeStuff = paste('TreeList_v',chromv,'_v',finev,sep="")
          new <- paste('new.split.matrix.',pre,'_v',chromv,'_v',finev,sep="")
          new.split.matrix <- get(new)
          maxsplit <- ncol(new.split.matrix)
          Tr <- get(treeStuff)
          if(!is.null(plotOrderRows)) {
            exclusion_samples <- rawdataRowOrder[!rawdataRowOrder%in%plotOrderRows]
          Tree =  drop.tip(Tr$ttree,tip=exclusion_samples)
          } else {
            Tree = Tr$ttree
          }
          tdend1 <- myapetodend(Tree,factor=1)
          if(!is.null(treeCut)) tdend <- cut(tdend1,h=treeCut)$upper else tdend=tdend1
                    
          treeOrderRows <- Tree$tip.label
          
          
      }

      if(tree==F){
          tdend=NULL
          if(!is.null(plotOrderRows)) treeOrderRows <- plotOrderRows[plotOrderRows%in%rawdataRowOrder] else treeOrderRows = rawdataRowOrder          
      }

        # get order for columns (rows in plot)
      if(!is.null(plotOrderCols)) treeOrderCols <- plotOrderCols[plotOrderCols%in%rawdataColOrder] else treeOrderCols = rawdataColOrder          
        
      print('getting row/column order')
          treeIndexRows <- match(treeOrderRows,SP_PR_HM_NA.attributes.all$Person.ID)
          treeDataIndexRows <- match(treeOrderRows,rawdataRowOrder,nomatch=0)
          treeIndexCols <- match(treeOrderCols,SP_PR_HM_NA.attributes.all$Person.ID)
          treeDataIndexCols <- match(treeOrderCols,rawdataColOrder,nomatch=0)

        # get labels for rows and columns
        print('getting labels')
                  Xlist <- niceLabels(labels[treeOrderRows]) # labels vector must have names 
                  Xtextcol <- "black"
                  Ylist <- niceLabels(labels[treeOrderCols]) 
                  Ytextcol <- "black"
                
        # get points/shapes (if labelPoints is specified, and matches labels)
          
        if(!is.null(labelPoints)){
          print('getting label points')
          
          labsX = unique(Xlist[Xlist!=""])
          labsY = unique(Ylist[Ylist!=""])
              
          if(length(which(!labsX%in%rownames(labelPoints)))>0){
            print("pointLabels don't have all the specified labels in it!")
            print(head(labsX))
            labelPoints=NULL
          } else {
              labPosX <- which(Xlist!="")
              labPosY <- which(Ylist!="")
    
              Xlist = paste(Xlist,"   ")
              Ylist = paste(Ylist,"   ")
                                    
              colsX=labelPoints[labsX,"colour"]
              shapesX=labelPoints[labsX,"shape"]
              colsY=labelPoints[labsY,"colour"]
              shapesY=labelPoints[labsY,"shape"]
          }
        }
              
        if(is.null(split_boundaries)) split_boundaries=c()
          
        
    # Start plotting
        for(matrixName in names(matrixTypes)){          
          
          print(paste('printing: ',matrixName))
          
              #get the matrix
              matrix = matrixTypes[matrixName]
              data1 <- as.matrix(get(paste(matrix,'_v',chromv,sep="")))
              
              #cap it if wanted
              if(cap==T) {
                print('capping matrix')
                cap=mean(data1,na.rm=T)+2*sd(data1,na.rm=T)
                data1[data1>cap] = cap
                print(cap)
              }              
              
              #make NAs 0
              data1[is.na(data1)] = 0
              print(range(data1))
          
          filename=paste(plotdir,'/',prefix,'.',matrixName,'_v',chromv,'_v',finev,'-%02d.png',sep="")
          
          png(filename,width=picSize[1],height=picSize[2],res=200)
                
          if(!is.null(labelPoints)) {             
            plotFinestructure(data1[treeDataIndexRows,treeDataIndexCols],labelsx=Xlist,labelsy=Ylist,dend=tdend,tickmarks=0,text.colx=Xtextcol,text.coly=Ytextcol,
                             cols=some.colorsEnd2,cex.axis=0.8,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8,t.col="transparent"),crt=45,
                             main=title) +
            abline(v=split_boundaries) + 
            points(rep(-10,length(labPosY)),labPosY,bg=colsY,pch=shapesY,col="black",xpd=NA)
            points(labPosX,rep(-10,length(labPosX)),bg=colsX,pch=shapesX,col="black",xpd=NA)
          } else {
            plotFinestructure(data1[treeDataIndexRows,treeDataIndexCols],labelsx=Xlist,labelsy=Ylist,dend=tdend,tickmarks=0,text.colx=Xtextcol,text.coly=Ytextcol,
                             cols=some.colorsEnd2,cex.axis=0.8,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8,t.col="transparent"),crt=45,
                             main=title) +
            abline(v=split_boundaries)
          }
                  #print(p)
        dev.off()
    }              
}

plotPaiwiseCoincidence <- function(matrix,chromv,plotdir,prefix,treeStuff=NULL,split.matrix=NULL,otherOrder=NULL,treeCut=7000,
                           plotOrderRows=NULL,plotOrderCols=NULL,
                           labelPoints=NULL,title="",labels=NULL,split_boundaries=NULL,cap=T,
                           some.colorsEnd2=NULL,picSize=c(2200,2000),filenameSuffix="",...){ 
  
  # get colours
    if(is.null(some.colorsEnd2)) some.colorsEnd2 <- colorRampPalette(c("yellow","red","blue","black"),interpolate="linear")(100) # as above, but with a dark grey final for capped values
  
  # get row/col order of raw matrices              
              data1 <- matrix             
              rawdataRowOrder = rownames(data1)
              rawdataColOrder = colnames(data1)

  # get order for rows (columns in plot)
      if(!is.null(treeStuff)){        
        print('getting tree info')
        # need to input Tree information for dendrogram (i.e the name of the TreeList file) + the new.split.matrix associated with the tree          
          new.split.matrix <- split.matrix
          maxsplit <- ncol(new.split.matrix)
          Tr <- get(treeStuff)
          if(!is.null(plotOrderRows)) {
            exclusion_samples <- Tr$ttree$tip.label[!Tr$ttree$tip.label%in%plotOrderRows]
            print(head(exclusion_samples))
          Tree =  drop.tip(Tr$ttree,tip=exclusion_samples)
          } else {
            Tree = Tr$ttree
          }
          tdend1 <- myapetodend(Tree,factor=1)
          if(!is.null(treeCut)) tdend <- cut(tdend1,h=treeCut)$upper else tdend=tdend1
                    
          treeOrderRows <- Tree$tip.label
          
          
      }

      if(is.null(treeStuff)){
          tdend=NULL
          if(!is.null(plotOrderRows)) treeOrderRows <- plotOrderRows[plotOrderRows%in%rawdataRowOrder] else treeOrderRows = rawdataRowOrder          
      }

        # get order for columns (rows in plot)
      if(!is.null(plotOrderCols)) treeOrderCols <- plotOrderCols[plotOrderCols%in%rawdataColOrder] else treeOrderCols = rawdataColOrder          
        
      print('getting row/column order')
          treeIndexRows <- match(treeOrderRows,SP_PR_HM_NA.attributes.all$Person.ID)
          treeDataIndexRows <- match(treeOrderRows,rawdataRowOrder,nomatch=0)
          treeIndexCols <- match(treeOrderCols,SP_PR_HM_NA.attributes.all$Person.ID)
          treeDataIndexCols <- match(treeOrderCols,rawdataColOrder,nomatch=0)

        # get labels for rows and columns
        print('getting labels')
        if(is.null(labels)) print('you need a named list of labels')
                  Xlist <- niceLabels(labels[treeOrderRows]) # labels vector must have names 
                  Xtextcol <- "black"
                  Ylist <- niceLabels(labels[treeOrderCols]) 
                  Ytextcol <- "black"
                
        # get points/shapes (if labelPoints is specified, and matches labels)
          
        if(!is.null(labelPoints)){
          print('getting label points')
          
          labsX = unique(Xlist[Xlist!=""])
          labsY = unique(Ylist[Ylist!=""])
              
          if(length(which(!labsX%in%rownames(labelPoints)))>0){
            print("pointLabels don't have all the specified labels in it!")
            print(head(labsX))
            labelPoints=NULL
          } else {
              labPosX <- which(Xlist!="")
              labPosY <- which(Ylist!="")
    
              Xlist = paste(Xlist,"   ")
              Ylist = paste(Ylist,"   ")
                                    
              colsX=labelPoints[labsX,"colour"]
              shapesX=labelPoints[labsX,"shape"]
              colsY=labelPoints[labsY,"colour"]
              shapesY=labelPoints[labsY,"shape"]
          }
        }
              
        if(is.null(split_boundaries)) split_boundaries=c()
          
        
    # Start plotting
 
              #cap it if wanted
              if(cap==T) {
                print('capping matrix')
                cap=mean(data1,na.rm=T)+2*sd(data1,na.rm=T)
                data1[data1>cap] = cap
                print(cap)
              }              
              
              #make NAs 0
              data1[is.na(data1)] = 0
              print(range(data1))
          
          filename=paste(plotdir,'/',prefix,'-',filenameSuffix,'.png',sep="")
          
          png(filename,width=picSize[1],height=picSize[2],res=200)
                
          if(!is.null(labelPoints)) {             
            plotFinestructure(data1[treeDataIndexRows,treeDataIndexCols],labelsx=Xlist,labelsy=Ylist,dend=tdend,tickmarks=0,text.colx=Xtextcol,text.coly=Ytextcol,
                             cols=some.colorsEnd2,cex.axis=0.8,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8,t.col="transparent"),crt=45,
                             main=title,...) +
            abline(v=split_boundaries) + 
            points(rep(-10,length(labPosY)),labPosY,bg=colsY,pch=shapesY,col="black",xpd=NA)
            points(labPosX,rep(-10,length(labPosX)),bg=colsX,pch=shapesX,col="black",xpd=NA)
          } else {
            plotFinestructure(data1[treeDataIndexRows,treeDataIndexCols],labelsx=Xlist,labelsy=Ylist,dend=tdend,tickmarks=0,text.colx=Xtextcol,text.coly=Ytextcol,
                             cols=some.colorsEnd2,cex.axis=0.8,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=0.8,t.col="transparent"),crt=45,
                             main=title,...) +
            abline(v=split_boundaries)
          }
                  #print(p)
        dev.off()
                       
}
  
  
  
###############
# Function to get date from generation time


gensToDates=function(gens,startDate,genTime=28){  
  startDate - genTime*(gens-1)  
}

################################## 
# Plotting functions for mixture model

barPlotSpreadOut <- function(data,subset="All",bootstraps=T,limit=0.001,...){

boot05 <- percentile05bootstrap[rownames(data),]
boot95 <- percentile95bootstrap[rownames(data),]
maxes <- apply(data,2,max)
maxes <- maxes[order(maxes,decreasing=T)]
these <- which(maxes>limit)
  
png(paste('context_plots/',version_string,"barPlots",subset,"-mixLimit",limit,".png",sep=""),height=1200,width=4000,res=200)
    #par(mfrow=c(1,1+length(these)))
    receivers <- c(1:nrow(data))
    if(bootstraps==T){
      xmax <- max(boot95[,names(these)])
      print(xmax)
    } else {
      xmax <- max(data[,names(these)])
    }

    #xmax=1
   # maxes <- apply(boot95[,these],2,max)
   # maxes <- c(3*max(maxes),maxes)
   # Widths <- maxes/sum(maxes)*100
    Widths <- c(3,rep(2,length(these)))
    
    layout(mat=matrix(c(1:(1+length(these))),nrow=1,ncol=(1+length(these))),widths=Widths)
    Labs <- sapply(rownames(data),getPopnLabel,popnLabels=popnLabels)[2,]
    par(mar=c(14,7,0,0),mgp=c(1,-4,0),cex=0.5)
    b <- barplot(rep(0,nrow(data)),horiz=T,col="transparent",xlim=c(0,1),las=2,names.arg=Labs,axes=F,border=F,xpd=NA)
    #dev.off()
    
    par(mar=c(14,0,0,0.3),mgp=c(3,0.8,-0.1),lwd=0.5,srt=300,...)
    for (i in names(these)){
      #xmax <- max(boot95[,i])
      donorLab <- getPopnLabel(i,popnLabels=popnLabels)[2]
      print(xmax)
    b <- barplot(data[,i],horiz=T,col=recieverLegend2[i,"colour"],xlim=c(0,xmax),las=2,names.arg=F,axes=F,xpd=NA)
    axis(1,cex=0.2,lwd=0.5)
      if(bootstraps==T){
           segments(boot05[receivers,i], b, boot95[receivers,i],  b,xpd=NA)
           segments(boot05[receivers,i], b - 0.1,  boot05[receivers,i],b + 0.1,xpd=NA)
           segments(boot95[receivers,i], b - 0.1,  boot95[receivers,i],b + 0.1,xpd=NA)
      }
    legend(x=0,y=-3,xpd=NA,legend=donorLab,col=recieverLegend2[i,"colour"],
           pch=as.numeric(recieverLegend2[i,"shape"]),
           pt.bg=recieverLegend2[i,"colour"],
           bty="n",
           horiz=T)
    
      }

dev.off()
}

barPlotMatrix <- function(data,subset="All",bootstraps=T,limit=0.001,dir='',
                          plotDonorLabels=T,height=NULL,width=NULL,Oma=c(12,8,1,0),
                          yTextMargin=2,xTextMargin=1,printLabs=T,labelStyle=2,
                          variableHeight=F,colourCols=TRUE,PDF=FALSE,diagText="x",...){

boot05 <- percentile05bootstrap[rownames(data),]
boot95 <- percentile95bootstrap[rownames(data),]
maxes <- apply(data,2,max)
#maxes <- maxes[order(maxes,decreasing=T)]
these <- which(maxes>limit)
  
#png(paste('context_plots',dir,'/',version_string,"barPlots",subset,"-mixLimit",limit,".png",sep=""),height=50*ncol(data),width=50*nrow(data),res=300)
if(is.null(height)) height = 50*length(these)+600
if(is.null(width)) width = 50*nrow(data)
  
if(!PDF) png(paste('context_plots',dir,'/',version_string,"barPlots",subset,"-mixLimit",limit,".png",sep=""),height=height,width=width,res=300)
if(PDF) pdf(paste('context_plots',dir,'/',version_string,"barPlots",subset,"-mixLimit",limit,".pdf",sep=""),height=height/300,width=width/300)

    #par(mfrow=c(1,1+length(these)))
    receivers <- c(1:nrow(data))
    if(bootstraps==T){
      xmax <- max(boot95[,names(these)])
      print(xmax)
    } else {
      xmax <- max(data[,names(these)])
    }

    #xmax=1
   # maxes <- apply(boot95[,these],2,max)
   # maxes <- c(3*max(maxes),maxes)
   # Widths <- maxes/sum(maxes)*100
        
    #layout(mat=matrix(c(1:(1+length(these))),ncol=1,nrow=(1+length(these))),heights=Heights)
    layout(mat=matrix(c(1:length(these)),ncol=1,nrow=length(these)))
    Labs <- sapply(rownames(data),getPopnLabel,popnLabels=popnLabels)[labelStyle,]
    #par(mar=c(14,7,0,0),mgp=c(1,-4,0),cex=0.5)
    #b <- barplot(rep(0,nrow(data)),horiz=F,col="transparent",xlim=c(0,1),las=2,names.arg=Labs,axes=F,border=F,xpd=NA)
    #dev.off()
  
    #par(mar=c(14,0,0,0.3),mgp=c(3,0.8,-0.1),lwd=0.5,srt=300,...)
    par(mar=c(0,0,0,0),mgp=c(0,0,0),lwd=0.2,srt=45,oma=c(10,7,1,0),cex=0.5)
   # par(...)
    par(oma=Oma,...)
    c = 1
    for (i in names(these)){
      #xmax <- max(boot95[,i])
      donorLab <- getPopnLabel(i,popnLabels=popnLabels)[labelStyle]
      print(xmax)
      if(colourCols) thisColour=recieverLegend2[rownames(data),"colour"] else thisColour = rep(recieverLegend2[i,"colour"],length(rownames(data)))
      doBorder=rep(NA,length(thisColour))
      doBorder[thisColour=="white"] = "black"
      if(variableHeight){
        xmax = max(data[,i])  # This actually sets the ymax limit
      }
      
      b <- barplot(data[,i],horiz=F,col=thisColour,ylim=c(0,xmax),las=2,names.arg=F,axes=F,xpd=NA,border=doBorder,width=1)
    #axis(2,cex=0.2,lwd=0.5)
      # Plot X on the diagonal
      text(x=b[colnames(data)==i],y=xmax/2,labels=diagText,srt=0,cex=0.8)
      
      if(bootstraps==T){
        nonZero = data[,i]!=0  # don't plot bootstraps for anything with a zero point estimate.
           segments(b[nonZero],boot05[receivers[nonZero],i],b[nonZero],boot95[receivers[nonZero],i],xpd=NA)
           segments(b[nonZero]-0.1,boot05[receivers[nonZero],i],b[nonZero]+0.1,boot05[receivers[nonZero],i],xpd=NA)
           segments(b[nonZero]-0.1,boot95[receivers[nonZero],i],b[nonZero]+0.1,boot95[receivers[nonZero],i],xpd=NA)
           if(sum(!nonZero)!=0) segments(b[!nonZero]-0.02,rep(0,sum(!nonZero)),b[!nonZero]+0.02,rep(0,sum(!nonZero)),xpd=NA) # a short horizontal line for zero.
      }
      
      print(xmax)
    a <- legend(x=-3*yTextMargin/2,y=xmax/2,xpd=NA,legend="",col="black",
           pch=as.numeric(recieverLegend2[i,"shape"]),
           pt.bg=recieverLegend2[i,"colour"],
           bty="n",
           horiz=T,cex=2,yjust=0.5)
      
    if(plotDonorLabels) text(a$text$x-yTextMargin,a$text$y,donorLab,xpd=NA,pos=2)
      
    # left-hand (surrogates) axis
    if(c!=1) axis(2,at = c(0,par("usr")[4]),labels=F,lwd=0.5,pos=c(0,par("usr")[1]),tcl=-0.3)    
    #if(c==1) axis(2,at = c(0,par("usr")[4]),lwd=0.8,line=-0.5,las=2)
    if(c==1) axis(2,at = c(0,par("usr")[4]),lwd=0.5,pos=c(0,par("usr")[1]),las=2,hadj=2,tcl=-0.3)
    
      c = c+1
      
      }
    text(b,-xTextMargin,Labs,xpd=NA,pos=2)
    points(b,rep(-xTextMargin/2,length(b)),
           pch=as.numeric(recieverLegend2[rownames(data),"shape"]),bg=recieverLegend2[rownames(data),"colour"],col="black",
           xpd=NA,cex=2)
    par(srt=0)

    
  if(printLabs){
    mtext("Ancestry profile components",side=2,xpd=NA,outer=T,line=5+yTextMargin,cex=0.8,col=add.alpha("black",0.5))
    #par(srt=0)
    mtext("Target groups",side=1,xpd=NA,outer=T,line=8+xTextMargin,cex=0.8,col=add.alpha("black",0.5))
  }
    par(opar)
dev.off()
}

#barPlotMatrix(data[donorOrder,donorOrder],subset="post-mix-Donors",bootstraps=T,limit=0.001,dir=directory,width=1500,height=1500,colourCols=TRUE,oma=c(11,9,1,0),cex=0.4,yTextMargin = 1.25)
#barPlotMatrix(data[donorOrder,donorOrder],subset="post-mix-Donors",bootstraps=T,limit=0.001,dir=directory,width=1500,height=1500,colourCols=TRUE,oma=c(11,9,1,0),cex=0.4,yTextMargin = 1.25,PDF=TRUE)

# boxplot of within-popn homogeneity

boxPlotByPopn = function(meanPopnCopydiffs,name="",dir = '',ylab="",filename=NULL,mar=c(15,5,2,2),...){
  
    popn = paste("donorPopn_",DONOR.split.matrix[names(meanPopnCopydiffs),2],sep="")
    popn[popn=="donorPopn_NA"] = paste("SpainPopn_",SPAIN.split.matrix[gsub("SPAIN_A2_","",names(meanPopnCopydiffs[popn=="donorPopn_NA"])),2],sep="")
    popn=factor(popn,levels=unique(popn))
    
    cols = recieverLegend2[levels(popn),"colour"]
    shapes = recieverLegend2[levels(popn),"shape"]
    
    if(is.null(filename)) filename=paste('context_plots',dir,'/',version_string,"boxPlot",name,sep="")
    
    #png(paste0(filename,".png"),height=800,width=1500,res=150)
    pdf(paste0(filename,".pdf"),height=800/150,width=1500/150,bg="transparent")
    
    par(mar=mar,...)    
      boxplot(meanPopnCopydiffs~popn,las=2,col=cols,outpch=shapes,outbg=cols,
            names=sapply(levels(popn),getPopnLabel,popnLabels)[1,],ylab=ylab)
    
    dev.off()
    par(opar)
}


piesOnMaps2=function(props,pieCentres,baseMap=NULL,pieColours=NULL,
                     radius=0.6,adj=c(1.3,1),dummyMap=commapborders,
                     mapLabels=commaplabels,moreText=NULL){
  
  #get base map and dummy map for adding pies to
  if(is.null(baseMap)) baseMap <- update(UnderMap2,main="") + update(commapborders2,col.regions=hsv(0,0,0.8),alpha.regions=0.5)
  p <- update(dummyMap,col.regions="transparent",col="transparent")
    
  #do you have centres?
  if(is.null(pieCentres)) {
    stop
    print('Where are the pies supposed to go, silly??')
  }
  pieNames=colnames(props)
     
  #get colours
  if(is.null(pieColours)) {
    get <- makeFactorLegend(rownames(props))
    pieColours <- get$colour
    names(pieColours) = get$Factor
  }
    
  #get pies
  pies <- lapply(pieNames,FUN=function(x) {  
            getPies2(props[,x],centre=pieCentres[x,c("X","Y")],polyID=x,radius=radius,adj=adj,colours=pieColours)
          })
  print(pies[[1]])
  
  #### put them on a map!
  pieMap = baseMap
  
  #map labels
  if(!is.null(mapLabels)) pieMap = pieMap +  update(p,sp.layout=list(mapLabels))
  
  #pies
  pieMap = pieMap + update(p,sp.layout=pies)
  
  #pie labels
  if(!is.null(moreText)) pieMap = pieMap + update(p,sp.layout=moreText)
    
  return(pieMap)
  
}

# get date labels for pies on maps. This requires pieCentres (the same as in piesOnMaps2), as well as Boots, and fitQuality objects
dateLabels <- function(thisNames,pieCentres,null=T,textcol="black"){
  toExtract = rownames(pieCentres[thisNames,])
  if("SpainPopn_49"%in%thisNames){
      toExtract[toExtract=="SpainPopn_49"]= "Portugal_6"
  }
  if(null){
    confDates <- sapply(paste(toExtract,".null",sep=""),function(x) quantile(Boots[[x]][,"date1.est.boot"],c(0.025,0.975)))
    singleDates <- fitQuality[paste("LOO-SP_PR_HM_NA.",paste(toExtract,".null",sep=""),sep=""),"gen.1date"]
  } else {
    confDates <- sapply(toExtract,function(x) quantile(Boots[[x]][,"date1.est.boot"],c(0.025,0.975)))
    singleDates <- fitQuality[paste("LOO-SP_PR_HM_NA.",toExtract,sep=""),"gen.1date"]   
  }
  
  confLabels <- list("sp.text",loc=cbind(pieCentres[thisNames,1],pieCentres[thisNames,2]-0.2),txt=paste(signif(startDate-29*confDates[2,],2),signif(startDate-29*confDates[1,],2),sep="-"),col=textcol,cex=0.5,alpha=0.9)
  singleLabels <- list("sp.text",loc=cbind(pieCentres[thisNames,1],pieCentres[thisNames,2]),txt=paste(signif(startDate-29*singleDates,2)),col=textcol,cex=1,alpha=0.6)

  return(list("conf"=confLabels,"single"=singleLabels))
}

############################## 
# Plotting functions for GT

PC.barplot.dates.Old <- function(PCinfo,title="",fullSurr=T,sources=NULL,rectanglesOnly=F,dateCol="gen.1date",bootstraps=F,bootsCol="date1.est.boot",
                             sampSizes=NULL,subsetPopns=NULL,years=F,genTime=29,plainNames=F,makeLegend=T,surrLeg=recieverLegend3,...){
          
          # get a and b-side + labels
          PCdataA <- PCinfo$a[,-1]
          PCdataA <- PCdataA[,gsub("LOO-SP_PR_HM_NA.","",rownames(fitQuality))]
          if(!is.null(subsetPopns)) PCdataA <- PCdataA[,colnames(PCdataA)%in%subsetPopns]
          #PCdataA <- PCdataA[,order(match(gsub('.null','',colnames(PCdataA)),ord))]
          xA <- apply(as.matrix(PCdataA),2,function(x) x[-length(x)]*x[length(x)])
          
          PCdataB <- PCinfo$b[,-1]
          PCdataB <- PCdataB[,colnames(PCdataA)]
          xB <- apply(as.matrix(PCdataB),2,function(x) x[-length(x)]*x[length(x)])
          
          dates <- fitQuality[match(paste("LOO-SP_PR_HM_NA.",colnames(PCdataA),sep=""),rownames(fitQuality)),]
          
          names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[2,]
          if(plainNames==T) names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[1,]
          names[grep("null",colnames(PCdataA))] <- paste(names[grep("null",colnames(PCdataA))],".null",sep="")
          leg <- surrLeg[gsub('.null','',colnames(PCdataA)),]
          
                    
          #layout(mat=matrix(c(1:length(colnames(PCdataA))),nrow=1,ncol=length(colnames(PCdataA))))
          par(...)
          
          barw=0.35
          barl=30
          linew= barw*1.3
          Ylabel="nGenerations"
          if(years==T)  Ylabel="Date (CE)"
            
          if(rectanglesOnly==F) {
            plot(c(1:ncol(PCdataA)),dates[,dateCol],ylim=c(min(dates[,dateCol],na.rm=T)-barl,max(dates[,dateCol],na.rm=T)+barl),bty="n",main=title,
               xaxt="n",yaxt="n",xlab="",col="transparent",ylab=Ylabel)
            if(years==T) axis(2,at=seq(0,max(dates[,dateCol],na.rm=T)+barl,10),labels=startDate - genTime*seq(0,max(dates[,dateCol],na.rm=T)+barl,10))
            if(years==F) axis(2,)
                    tck <- axis(1,at=c(1:ncol(PCdataA)),labels=F,tick=F)
                    text(tck,par("usr")[3]-2,labels=names,srt=60,pos=2,xpd=NA,cex=0.8,offset=c(-0.3,0))
          }
          #b <- barplot(dates$gen.1date,col="gray",border=NA,width=10,names.arg=names,las=2)
          
          for (pop in c(1:ncol(PCdataA))){
                
                date <- dates[pop,dateCol]
                PCa <- xA[,pop]
                PCb <- xB[,pop]
                
                if(fullSurr==F){
                  PCa <- PCdataA[nrow(PCdataA),pop]
                  PCb <- PCdataB[nrow(PCdataB),pop]
                  names(PCa) <- dates[pop,sources[1]] 
                  names(PCb) <- dates[pop,sources[2]]
                                            
                }
                if(rectanglesOnly==F) rect(pop-barw,par("usr")[3],pop+barw,par("usr")[4],col="lightgray",border=NA)
                
                #A
                this=0
                for (sur in names(PCa)[PCa>0]){
                        
                        rect(pop-barw,date+this,pop+barw,date+this+PCa[sur]*barl,col=surrLeg[sur,"colour"])
                        this <- PCa[sur]*barl + this      
                }
                #B
                this=0
                for (sur in names(PCb)[PCb>0]){
                        rect(pop-barw,date-this,pop+barw,date-this-PCb[sur]*barl,col=surrLeg[sur,"colour"])
                        this <- PCb[sur]*barl + this
                }
                
                lines(c(pop-linew,pop+linew),c(date,date),lwd=2)
                
                if(bootstraps==T){
                  booties <- quantile(Boots[[colnames(PCdataA)[pop]]][,bootsCol],c(0.05,0.95))
                  
                  lines(c(pop,pop),c(booties[1],booties[2]),lwd=2)
                  lines(c(pop-linew,pop+linew),c(booties[1],booties[1]),lwd=2)
                  lines(c(pop-linew,pop+linew),c(booties[2],booties[2]),lwd=2)
                }
          
          }
          
          if(rectanglesOnly==F){ 
            
              points(tck,rep(par("usr")[3]+2,ncol(PCdataA)),pch=leg$shape,bg=leg$colour,xpd=NA)
              if(!is.null(sampSizes)){
                sizeLabs <- sampSizes[gsub(".null","",colnames(PCdataA))]
                text(tck,rep(par("usr")[3]+4,ncol(PCdataA)),labels=sizeLabs,xpd=NA,cex=0.3)
              }
              nonzerosurr <- which(surrLeg$Factor%in%unique(c(rownames(PCdataA)[rowSums(PCdataA)>0.01],rownames(PCdataB)[rowSums(PCdataB)>0])))
              if(fullSurr==F){
                nonzerosurr <- surrLeg$Factor%in%unique(c(dates[,sources[1]],dates[,sources[2]]))
              }
              surrNames <- sapply(surrLeg$Factor,getPopnLabel,popnLabels)[2,]
          }
        if(makeLegend==T) legend("topright",legend=surrNames[nonzerosurr],col=surrLeg$colour[nonzerosurr],fill=surrLeg$colour[nonzerosurr],border=NA,bg="white",cex=0.8,inset=c(-0.1,0.05),xpd=NA)
  }



PC.barplot.dates <- function(PCinfo,title="",fullSurr=T,sources=NULL,rectanglesOnly=F,dateCol=dateCol,bootstraps=F,bootsCol="date1.est.boot",
                             sampSizes=NULL,conclusions=NULL,subsetPopns=NULL,years=F,genTime=29,popOrder=NULL,plainNames=F,makeLegend=T,
                             surrLeg=recieverLegend3,lines=NULL,fixYaxis=NULL,legMax=0.02,
                             ...){
          
          # get a and b-side + labels
          PCdataA <- PCinfo$a[,-1]
          PCdataA <- PCdataA[,gsub("LOO-SP_PR_HM_NA.","",rownames(fitQuality))]
          if(!is.null(subsetPopns)) PCdataA <- PCdataA[,colnames(PCdataA)%in%subsetPopns]
          if(!is.null(popOrder)) PCdataA <- PCdataA[,popOrder]
          
          #PCdataA <- PCdataA[,order(match(gsub('.null','',colnames(PCdataA)),ord))]
          xA <- apply(as.matrix(PCdataA),2,function(x) x[-length(x)]*x[length(x)])
          
          PCdataB <- PCinfo$b[,-1]
          PCdataB <- PCdataB[,colnames(PCdataA)]
          xB <- apply(as.matrix(PCdataB),2,function(x) x[-length(x)]*x[length(x)])
          
          dates <- fitQuality[match(paste("LOO-SP_PR_HM_NA.",colnames(PCdataA),sep=""),rownames(fitQuality)),]
          
          names <- sapply(gsub('.null|.2dates','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[3,]
          if(plainNames==T) names <- sapply(gsub('.null|.2dates','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[1,]
          names[grep("null",colnames(PCdataA))] <- paste(names[grep("null",colnames(PCdataA))],".null",sep="")
          leg <- surrLeg[gsub('.null','',names),]
          
                    
          #layout(mat=matrix(c(1:length(colnames(PCdataA))),nrow=1,ncol=length(colnames(PCdataA))))
          par(...)
          
          barw=0.35
          barl=30
          linew= barw*1.5
          Ylabel="nGenerations"
          if(years==T)  Ylabel="Date (CE)"
            
          if(rectanglesOnly==F) {
            #plot(c(1:ncol(PCdataA)),dates[,dateCol],ylim=c(min(dates[,dateCol],na.rm=T)-barl,max(dates[,dateCol],na.rm=T)+barl),bty="n",main=title,
             #  xaxt="n",yaxt="n",xlab="",col="transparent",ylab=Ylabel)
            if(is.null(fixYaxis)) ymax=max(dates[,dateCol],na.rm=T)+barl else ymax=fixYaxis
            
            plot(c(1:ncol(PCdataA)),dates[,dateCol],ylim=c(-10,ymax),bty="n",main=title,
               xaxt="n",yaxt="n",xlab="",col="transparent",ylab=Ylabel)
            if(years==T) axis(2,at=seq(0,max(dates[,dateCol],na.rm=T)+barl,10),labels=startDate - genTime*seq(0,max(dates[,dateCol],na.rm=T)+barl,10))
            if(years==F) axis(2,)
                    tck <- axis(1,at=c(1:ncol(PCdataA)),labels=F,tick=F)
                    text(tck,par("usr")[3]-2,labels=names,srt=60,pos=2,xpd=NA,cex=0.8,offset=c(-0.3,0))
          }
          #b <- barplot(dates$gen.1date,col="gray",border=NA,width=10,names.arg=names,las=2)
          
          barl = par("usr")[4] - par("usr")[3]
          
          #Order surrogate proportions by average size
          forOrd=rowMeans(xA,na.rm=T)
          surrOrderA = names(forOrd)[forOrd>0]
          surrOrderA = surrOrderA[order(forOrd[forOrd>0])]
          
          forOrd=rowMeans(xB,na.rm=T)
          surrOrderB = names(forOrd)[forOrd>0]
          surrOrderB = surrOrderB[order(forOrd[forOrd>0])]
          
          #surrLeg$gtColours = NA
          #surrLeg[surrOrderB,"gtColours"] = colorRampPalette(c("green","blue","purple"))(length(surrOrderB))
          #surrLeg[surrOrderA,"gtColours"] = colorRampPalette(c("white","yellow","red"))(length(surrOrderA))
          
          for (pop in c(1:ncol(PCdataA))){
                
                date <- dates[pop,dateCol]
                PCa <- xA[,pop]
                PCb <- xB[,pop]
                
                if(fullSurr==F){
                  PCa <- PCdataA[nrow(PCdataA),pop]
                  PCb <- PCdataB[nrow(PCdataB),pop]
                  names(PCa) <- dates[pop,sources[1]] 
                  names(PCb) <- dates[pop,sources[2]]
                  surrOrderA = names(PCa)     
                  surrOrderB = names(PCb)
                }
                if(rectanglesOnly==F) rect(pop-barw,par("usr")[3],pop+barw,par("usr")[4],col="lightgray",border=NA)
                
                #A                
                par(xpd=NA)
                this=0
                gap=0.03
                newBarl=par("usr")[4]*(1-gap)

                for (sur in surrOrderA){
                        
                        rect(pop-barw,this,pop+barw,this+PCa[sur]*newBarl,col=surrLeg[sur,"colour"])
                        this <- PCa[sur]*newBarl + this      
                }
                #B
                this=this + gap*par("usr")[4]
                
                for (sur in surrOrderB){
                  
                        rect(pop-barw,this,pop+barw,this+PCb[sur]*newBarl,col=surrLeg[sur,"colour"])
                        this <- PCb[sur]*newBarl + this
                }
                
                #print(c(pop-linew,pop+linew))
                #print(c(date,date))
                
                graphics::lines(c(pop-linew,pop+linew),c(date,date),lwd=3)
                
                if(bootstraps==T){
                  print('plotting bootstraps')
                  if(dateCol%in%c("gen.2dates.date1","gen.2dates.date2")) {
                          booties <- quantile(Boots[[paste0(colnames(PCdataA)[pop],".2dates")]][,bootsCol],c(0.05,0.95))                  
                    } else {
                          booties <- quantile(Boots[[colnames(PCdataA)[pop]]][,bootsCol],c(0.05,0.95))            
                    }
                  graphics::lines(c(pop,pop),c(booties[1],booties[2]),lwd=2)
                  graphics::lines(c(pop-linew,pop+linew),c(booties[1],booties[1]),lwd=2)
                  graphics::lines(c(pop-linew,pop+linew),c(booties[2],booties[2]),lwd=2)
                }
          
          }
          
          if(rectanglesOnly==F){ 
            
              points(tck,rep(par("usr")[3]+2,ncol(PCdataA)),pch=leg$shape,bg=leg$colour,xpd=NA)
              
              #add conclusions legend
              if(!is.null(conclusions)){   
                print('printing conclusions')
                concl=conclusions[colnames(PCdataA)]              
                concLeg=c(4,8,7)
                names(concLeg)=c("one-date","one-date-multiway","multiple-dates")
                print(concLeg[concl])
                points(tck,rep(par("usr")[3]+7,ncol(PCdataA)),pch=concLeg[concl],col="black",xpd=NA,lwd=1)
                legend("bottomright",legend=unique(concl),pch=concLeg[unique(concl)],col="black",inset=c(-0.25,0))
              }
              
              if(!is.null(sampSizes)){
                sizeLabs <- sampSizes[gsub(".null","",colnames(PCdataA))]
                text(tck,rep(par("usr")[3]+4,ncol(PCdataA)),labels=sizeLabs,xpd=NA,cex=0.3)
              }
              maxesA=apply(PCdataA,1,max,na.rm=T)
              maxesB=apply(PCdataB,1,max,na.rm=T)
              nonzerosurr <- which(surrLeg$Factor%in%unique(c(names(maxesA)[maxesA*maxesA["totalPCprop"]>legMax],names(maxesB)[maxesB*maxesB["totalPCprop"]>0.05])))
              if(fullSurr==F){
                nonzerosurr <- surrLeg$Factor%in%unique(c(dates[,sources[1]],dates[,sources[2]]))
              }
              
              
              #surrNames <- sapply(surrLeg$Factor,getPopnLabel,popnLabels)[2,]
          }
          
        if(makeLegend==T) legend("topright",legend=surrLeg$Factor[nonzerosurr],col=surrLeg$colour[nonzerosurr],fill=surrLeg$colour[nonzerosurr],border=NA,bg="white",cex=0.8,inset=c(-0.25,0.05),xpd=NA)
        
        if(!is.null(lines)){
          
          par(xpd=F)
          abline(h=lines,lty=2,lwd=2,col="gray")
        }
    }


# testing above
#PCinfo <- PCset[[type]]
#filename=paste(plotdir,'/GTbarplots-',version,'-',type,'-%02d--test.png',sep="")
#popOrder=grep("null",colnames(PCinfo[[1]]),invert=T,value=T)
#popOrder = c("Portugal_6",popOrder[-28])
#popOrder = SpainGeogOrder[c(12,1:11,13,14,20,19,18,23:25,26:28,16,22,15,17,21)]
##popOrder = rownames(mixmat)[grep("Spain",rownames(mixmat))]

#png(filename,width=3000,height=2000,res=200)
#  par(opar)

#PC.barplot.dates(PCinfo=PCinfo,subsetPopns=popOrder,popOrder=popOrder,lines=c(24,42),title="1 date, 1st event, source mixture",sampSizes=sizes,
#                   mar=c(12,5,3,15),cex.axis=0.5,lwd=0.5,cex=1,bootstraps=F,bootsCol="date1.est.boot")
#dev.off()
# testing above

## new version of GT plots

plotGT2Old <- function(PCinfo,title="",dateCol="gen.1date",bootstraps=F,bootsCol="date1.est.boot",
                             sampSizes=NULL,conclusions=NULL,subsetPopns=NULL,years=F,genTime=29,popOrder=NULL,plainNames=F,makeLegend=T,
                             surrLeg=recieverLegend3,lines=NULL,fixXaxis=NULL,tol=0.01,fade=0.6,extraLines=NULL,dateLineWidth=8,
                             ...){
                    
          PCdataA <- PCinfo$a[,-1]
          PCdataA <- PCdataA[,gsub("LOO-SP_PR_HM_NA.","",rownames(fitQuality))]
          if(!is.null(subsetPopns)) PCdataA <- PCdataA[,colnames(PCdataA)%in%subsetPopns]
          if(!is.null(popOrder)) PCdataA <- PCdataA[,popOrder]
          
          #PCdataA <- PCdataA[,order(match(gsub('.null','',colnames(PCdataA)),ord))]
          xA <- apply(as.matrix(PCdataA),2,function(x) x[-length(x)]*x[length(x)])
          
          PCdataB <- PCinfo$b[,-1]
          PCdataB <- PCdataB[,colnames(PCdataA)]
          xB <- apply(as.matrix(PCdataB),2,function(x) x[-length(x)]*x[length(x)])
          
          dates <- fitQuality[match(paste("LOO-SP_PR_HM_NA.",colnames(PCdataA),sep=""),rownames(fitQuality)),]
          
          names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[2,]
          if(plainNames==T) names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[1,]
          names[grep("null",colnames(PCdataA))] <- paste(names[grep("null",colnames(PCdataA))],".null",sep="")
          leg <- surrLeg[gsub('.null','',names),]
          
                 
          if(is.null(fixXaxis)) xmax=max(dates[,dateCol],na.rm=T) else xmax=fixXaxis
          if(bootstraps){
                                                 
                  booties <- sapply(1:ncol(PCdataA),function(x) quantile(Boots[[colnames(PCdataA)[x]]][,bootsCol],c(0.025,0.975)))
                                    
                  if(is.null(fixXaxis)) {
                    xmin=min(booties,na.rm=T) 
                    xmax=max(booties,na.rm=T) 
                  }  else {
                    xmin=fixXaxis[1]
                    xmax=fixXaxis[2]
                  }
           }
          
  nTargets = dim(PCdataA)[2]
  
  
  par(bg="transparent",...)
            
  # plot surrogate components
  allDonorsA = names(which(rowSums((xA>tol),na.rm=T)>0))       
  allDonorsB = names(which(rowSums((xB>tol),na.rm=T)>0)) 
          
  donorOrderA = allDonorsA[order(rowSums(xA[allDonorsA,],na.rm=T))]
  donorOrderB = rev(allDonorsB[order(rowSums(xB[allDonorsB,],na.rm=T))])
          
          
  #layout(matrix(c(2,3,1), 1, 3, byrow=TRUE),widths=c((1-dateWidth)/2,(1-dateWidth)/2),dateWidth)
  par(mar=c(5,0,1,0))
  
  plot(NULL,xlim=c(xmin,xmax),ylim=c(1,nTargets)*dateLineWidth,yaxt="n",ylab=NA,
       bty="n",xlab="number of generations before present")
  
  
  xMin = par("usr")[1]
  xMax = par("usr")[2] 
  
          
  gap = (xMax-xMin)/30
          #vectors for widths of notches and colours
  yPos0 = (1:nTargets)*dateLineWidth-dateLineWidth/2+dateLineWidth/8
  yPos1 = (1:nTargets)*dateLineWidth+dateLineWidth/2-dateLineWidth/8
          
  yMax = max(yPos1)
          
   # side A
  isOverTol = xA[donorOrderA,]>tol   
  rectCols = isOverTol
  for(d in donorOrderA){
      rectCols[d,] = add.alpha(surrLeg[d,"colour"],fade)
  }
  rectCols[!isOverTol] = "lightgray"               
  #plot(NULL,xlim=c(1,length(donorOrderA)+1),ylim=c(1,nTargets)*dateLineWidth,yaxt="n",xaxt="n",xlab=NA,bty="n",ylab="minor side")          
  this = xmin
  step = ((xmax-xmin)-gap)/(length(donorOrderA)+length(donorOrderB))          
  for(donor in donorOrderA){
      rect(rep(this,nTargets),
                     yPos0,
                     rep(this+step,nTargets),
                     yPos1,
                     col=rectCols[donor,],border=NA)
      this = this+step
  }
          
          # side B
  isOverTol = xB[donorOrderB,]>tol   
  rectCols = isOverTol
  for(d in donorOrderB){
      rectCols[d,] = My.add.alpha(surrLeg[d,"colour"],fade)
  }
  rectCols[!isOverTol] = "lightgray"   
                      
  #plot(NULL,xlim=c(1,length(donorOrderB)+1),ylim=c(1,nTargets)*dateLineWidth,yaxt="n",xaxt="n",xlab=NA,bty="n",ylab="minor side")          
  this = this + gap
  for(donor in donorOrderB){
      graphics::rect(rep(this,nTargets),
                     yPos0,
                     rep(this+step,nTargets),
                     yPos1,
                     col=rectCols[donor,],border=NA)
      this = this+step
  }
                 
  
  # plot dates          
  dateWidth=0.8
  
            graphics::segments(dates[,dateCol],
                     yPos0,
                     dates[,dateCol],
                     yPos1,
                     lwd=3)
  
          
  if(bootstraps){
    print('plotting bootstraps') 
    graphics::segments(booties[1,],(1:nTargets)*dateLineWidth,booties[2,],(1:nTargets)*dateLineWidth,lwd=2)
  
  }   
  
  
          
          # labels
          points(x=xmin + (1:length(donorOrderA)*step)-step/2,y=rep(yMax+yMax/30,length(donorOrderA)),
           pch=surrLeg[donorOrderA,"shape"],col="black",bg=surrLeg[donorOrderA,"colour"],xpd=NA,cex=1.5)
          points(x=xmin + (length(donorOrderA)*step) +gap + (1:length(donorOrderB)*step)-step/2,y=rep(yMax+yMax/30,length(donorOrderB)),
           pch=surrLeg[donorOrderB,"shape"],col="black",bg=surrLeg[donorOrderB,"colour"],xpd=NA,cex=1.5)
             
          #text(x=60,y=yPos1,labels=colnames(rectCols))
          
             points(x=rep(xmax+xmax/50,nTargets),y=yPos1,pch=surrLeg[gsub(".null","",colnames(rectCols)),"shape"],
                bg=surrLeg[gsub(".null","",colnames(rectCols)),"colour"])                          
          
}



 
plotGT2 <- function(PCinfo,title="",dateCol="gen.1date",bootstraps=F,bootsCol="date1.est.boot",
                             sampSizes=NULL,conclusions=NULL,subsetPopns=NULL,years=F,rat=0.3,genTime=29,popOrder=NULL,plainNames=F,makeLegend=T,
                             surrLeg=recieverLegend3,bgcols=NULL,lines=NULL,fixXaxis=NULL,tol=0.01,applyFading=TRUE,propText=TRUE,fade=1,extraLines=NULL,dateLineWidth=8,textLabels=T,pointLabels=T,
                             mixmat=NULL,boot95=NULL,boot05=NULL,mixmatDonor="donorPopn_1",titleLine=1,titleMargin=5,
                             ...){
  
          PCdataA <- PCinfo$a[,-1]
          PCdataA <- PCdataA[,gsub("LOO-SP_PR_HM_NA.","",rownames(fitQuality))]
          if(!is.null(subsetPopns)) PCdataA <- PCdataA[,colnames(PCdataA)%in%subsetPopns]
          if(!is.null(popOrder)) PCdataA <- PCdataA[,popOrder]
          
          #PCdataA <- PCdataA[,order(match(gsub('.null','',colnames(PCdataA)),ord))]
          xA <- apply(as.matrix(PCdataA),2,function(x) x[-length(x)]*x[length(x)])
          
          PCdataB <- PCinfo$b[,-1]
          PCdataB <- PCdataB[,colnames(PCdataA)]
          xB <- apply(as.matrix(PCdataB),2,function(x) x[-length(x)]*x[length(x)])
          
          dates <- fitQuality[match(paste("LOO-SP_PR_HM_NA.",colnames(PCdataA),sep=""),rownames(fitQuality)),]
          
          names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[3,]
          if(plainNames==T) names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[1,]
          names[grep("null",colnames(PCdataA))] <- paste(names[grep("null",colnames(PCdataA))],".null",sep="")
          leg <- surrLeg[gsub('.null','',names),]
          
                 
          if(is.null(fixXaxis)) xmax=max(dates[,dateCol],na.rm=T) else xmax=fixXaxis[2]
          xmin=0
          if(bootstraps){
                                                 
                 if(dateCol %in% c("gen.2dates.date1","gen.2dates.date2")) {
                   booties <- sapply(1:ncol(PCdataA),function(x) quantile(Boots[[paste0(colnames(PCdataA)[x],".2dates")]][,bootsCol],c(0.025,0.975)))
                 } else {
                   booties <- sapply(1:ncol(PCdataA),function(x) quantile(Boots[[colnames(PCdataA)[x]]][,bootsCol],c(0.025,0.975)))
                 }
                  print(booties)
                  
                  if(is.null(fixXaxis)) {
                    xmin=min(booties,na.rm=T) 
                    xmax=max(booties,na.rm=T) 
                  }  else {
                    xmin=fixXaxis[1]
                    xmax=fixXaxis[2]
                  }
           }
          
  nTargets = dim(PCdataA)[2]
  popOrderNoNull = gsub(".null","",popOrder)
        
  par(bg="transparent",...)
            
  print(par()$mgp)
  #allDonorsA = names(which(rowSums((xA>tol),na.rm=T)>0))       
  #allDonorsB = names(which(rowSums((xB>tol),na.rm=T)>0)) 
  #allDonorsA = names(which(rowSums((xA>0),na.rm=T)>0))       
  #allDonorsB = names(which(rowSums((xB>0),na.rm=T)>0)) 
  allDonorsA = names(which(apply(PCdataA[-nrow(PCdataA),],1,function(x) sum(x>tol)>0)))       # at least one target above tolerance
  allDonorsB = names(which(apply(PCdataB[-nrow(PCdataB),],1,function(x) sum(x>tol)>0))) 
          
  donorOrderA = allDonorsA[order(rowSums(xA[allDonorsA,],na.rm=T))]
  donorOrderB = rev(allDonorsB[order(rowSums(xB[allDonorsB,],na.rm=T))])
          
          #vectors for widths of notches and colours
        yPos0 = (1:nTargets)*dateLineWidth-dateLineWidth/2+dateLineWidth/8
        yPos1 = (1:nTargets)*dateLineWidth+dateLineWidth/2-dateLineWidth/8
                
        yMax = max(yPos1)
        
  # set layout
  layout(t(c(1,2)),widths=c(rat,1-rat))
        
        
  if(!is.null(mixmat)){
        # plot mixmat components (made in mixtureModel with coancestry)
        mixMatData = mixmat[popOrderNoNull,mixmatDonor]
    
        par(mar=c(3,2,titleMargin,1))
        
        #plot(NULL,xlim=c(0,10),ylim=c(0,nTargets)*dateLineWidth+dateLineWidth/2-dateLineWidth/8,yaxt="n",ylab="Admixing groups",
        #     bty="n",xaxt="n")
        
        #xMin = par("usr")[1]
        #xMax = par("usr")[2] 
         
        b = barplot(mixMatData,horiz=TRUE,col=surrLeg[mixmatDonor,"colour"],
                    names=NA,
                    xlim=c(xmin,xmax),xpd=NA)
        print(c(xmin,xmax))
        title(main=paste0(getPopnLabel(mixmatDonor,popnLabels)[1]," fraction in ancestry profiles"),line=titleLine)
      
        segments(boot05[popOrderNoNull,mixmatDonor],b,boot95[popOrderNoNull,mixmatDonor],b,xpd=NA)
        segments(boot05[popOrderNoNull,mixmatDonor],b-0.1,boot05[popOrderNoNull,mixmatDonor],b+0.1,xpd=NA)
        segments(boot95[popOrderNoNull,mixmatDonor],b-0.1,boot95[popOrderNoNull,mixmatDonor],b+0.1,xpd=NA)            
  
        # point labels for target groups
        bgCols = surrLeg[popOrderNoNull,"colour"]
        pointBorderCols=rep("black",length(bgCols))
        
        points(x=rep(-1/50,nTargets),y=b,pch=surrLeg[popOrderNoNull,"shape"],xpd=NA,
              bg=bgCols,col=pointBorderCols)
          
        
    } else {
        
        # plot surrogate components
        par(mar=c(3,2,titleMargin,0))
        
        plot(NULL,xlim=c(0,10),ylim=c(0,nTargets)*dateLineWidth+dateLineWidth/2-dateLineWidth/8,yaxt="n",ylab="Admixing groups",
             bty="n",xaxt="n")
        
        xMin = par("usr")[1]
        xMax = par("usr")[2] 
                  
        gap = 0.5
        #if(propText) gap = 1
        
         # side A
        #isOverTol = xA[donorOrderA,]>tol   
        isOverTol = PCdataA[donorOrderA,]>tol   
        rectCols = isOverTol
        for(d in donorOrderA){
          d2 = paste0("donorPopn_",gsub("[^[:digit:]]","",d))
          print(d2)
          rectCols[d,] = add.alpha(surrLeg[d2,"colour"],fade)
        }
              
        rectCols[!isOverTol] = "lightgray"
        
        if(applyFading){
          rectCols = t(sapply(donorOrderA,FUN=function(i) My.add.alpha2(surrLeg[paste0("donorPopn_",gsub("[^[:digit:]]","",i)),"colour"],alpha=PCdataA[i,])))
          colnames(rectCols) = colnames(isOverTol) 
          rownames(rectCols) = rownames(isOverTol) 
        }
          
                
        this = 0
        step = (10-gap)/(length(donorOrderA)+length(donorOrderB))          
        for(donor in donorOrderA){
            graphics::rect(rep(this,nTargets),
                           yPos0,
                           rep(this+step,nTargets),
                           yPos1,
                           col=rectCols[donor,],border=NA)
            this = this+step    
        }
        # add text for total amount
        #middle = this/2
        middle = 0-0.1*this
        if(propText) text(x=middle,y=yPos0+(yPos1[1]-yPos0[1])/2,labels=PCdataA["totalPCprop",], xpd=NA)    
        if(!propText) mtext(text=paste0(paste0(100*range(PCdataA["totalPCprop",]),collapse=" - "),"%"),side=1,at=middle)          
                
        # draw line between side A and B
        middleLine = this + gap/2
        abline(v=middleLine,lty=3,lwd=2)
        
                # side B
        #isOverTol = xB[donorOrderB,]>tol   
        isOverTol = PCdataB[donorOrderB,]>tol     
        rectCols = isOverTol
        for(d in donorOrderB){
          d2 = paste0("donorPopn_",gsub("[^[:digit:]]","",d))
            rectCols[d,] = My.add.alpha(surrLeg[d,"colour"],fade)
        }
        rectCols[!isOverTol] = "lightgray"   
                            
        if(applyFading){
          rectCols = t(sapply(donorOrderB,FUN=function(i) My.add.alpha2(surrLeg[paste0("donorPopn_",gsub("[^[:digit:]]","",i)),"colour"],alpha=PCdataB[i,])))
          colnames(rectCols) = colnames(isOverTol) 
          rownames(rectCols) = rownames(isOverTol)   
        }
          
        #plot(NULL,xlim=c(1,length(donorOrderB)+1),ylim=c(1,nTargets)*dateLineWidth,yaxt="n",xaxt="n",xlab=NA,bty="n",ylab="minor side")          
        this0 = this + gap
        this=this0
        for(donor in donorOrderB){
            graphics::rect(rep(this,nTargets),
                           yPos0,
                           rep(this+step,nTargets),
                           yPos1,
                           col=rectCols[donor,],border=NA)
            this = this+step
        }
        # add text for total amount
        # if(propText) text(x=this - (this-this0)/2,y=yPos0+(yPos1[1]-yPos0[1])/2,labels=PCdataB["totalPCprop",])  
        #mtext(text=paste0(paste0(100*range(PCdataB["totalPCprop",]),collapse=" - "),"%"),side=1,at=this - middle)          
        
        # labels for donor populations
        if(textLabels){        
          # do text labels
          labsA = donorOrderA
          labsB = donorOrderB
                if(plainNames) {
                  #labsA = gsub("_|[[:digit:]]","",labsA); labsA = gsub("M-A-L","Mixed",labsA)
                  #labsB = gsub("_|[[:digit:]]","",labsB); labsB = gsub("M-A-L","Mixed",labsB)
                  labsA = sapply(paste0("donorPopn_",gsub("[^[:digit:]]","",labsA)),getPopnLabel,popnLabels)[1,]
                  labsB = sapply(paste0("donorPopn_",gsub("[^[:digit:]]","",labsB)),getPopnLabel,popnLabels)[1,]
            
                }
                if(pointLabels) fraqUp = 15 else fraqUp = 30
          
                text(x=seq(step,length(donorOrderA)*step,step)-step/2,y=rep(yMax+yMax/fraqUp,length(donorOrderA)),offset=step/4,
                     labels = labsA,xpd=NA,srt=45,pos=4,cex=0.7)
                text(x=length(donorOrderA)*step +gap + seq(step,length(donorOrderB)*step,step)-step/2,y=rep(yMax+yMax/fraqUp,length(donorOrderB)),offset=step/4,
                     labels = labsB,xpd=NA,srt=45,pos=4,cex=0.7)
        }
        if(pointLabels)  {
          # draw points labels
          points(x=seq(step,length(donorOrderA)*step,step)-step/2,y=rep(yMax+yMax/30,length(donorOrderA)),
                 pch=surrLeg[paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderA)),"shape"],col="black",bg=surrLeg[paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderA)),"colour"],xpd=NA,cex=1.5)
                points(x=length(donorOrderA)*step +gap + seq(step,length(paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderB)))*step,step)-step/2,y=rep(yMax+yMax/30,length(donorOrderB)),
                 pch=surrLeg[paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderB)),"shape"],col="black",bg=surrLeg[paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderB)),"colour"],xpd=NA,cex=1.5)          
        }               
        
                
    }
    
  # plot dates
  par(mar=c(3,0,titleMargin,6))
          
  print(par()$mgp)
  print( c(xmin,xmax) )
  
  plot(NULL,xlim=c(xmin,xmax),ylim=c(0,nTargets)*dateLineWidth+dateLineWidth/2-dateLineWidth/8,yaxt="n",ylab=NA,
       bty="n",xaxt="n")
  title(main="Admixture dates",line=titleLine)
  n = length(dateCol)      
  if(is.null(bgcols)) bgcols = surrLeg[popOrderNoNull,"colour"]  
  rect(rep(rep(xmin,n),n), yPos0, rep(xmax,n), yPos1,border=NA, col = bgcols)
  
  if(years) axis(1,at=seq(0,round(xmax,-1),10),labels=startDate - genTime*(seq(0,round(xmax,-1),10)-1))
  if(!years) axis(1)
  
  bootstrapLineCols = rep("black",nTargets)
  bootstrapLineCols[bgcols=="black"]="lightgray"
    
  dateWidth=0.8
  graphics::segments(dates[,dateCol],
                     yPos0,
                     dates[,dateCol],
                     yPos1,
                     lwd=3,col=bootstrapLineCols)
  
          
  if(bootstraps){
    print('plotting bootstraps') 
    graphics::segments(booties[1,],(1:nTargets)*dateLineWidth,booties[2,],(1:nTargets)*dateLineWidth,lwd=2,col=bootstrapLineCols)
  
  }   
  # point labels for target groups
  bgCols = surrLeg[popOrderNoNull,"colour"]
  pointBorderCols=rep("black",length(bgCols))
  dark = rgb2hsv(col2rgb(bgCols))["v",]<0.56
  pointBorderCols[dark] = "gray"

    points(x=rep(xmax-xmax/50,nTargets),y=yPos1-(yPos1[1]-yPos0[1])/2,pch=surrLeg[popOrderNoNull,"shape"],
        bg=bgCols,col=pointBorderCols)   
    
    # word labels for target groups
    print(popOrderNoNull)
    treeLabels = sapply(popOrderNoNull,getPopnLabel,popnLabels)[1,]
    text(x= rep(xmax+xmax/70,nTargets),y=yPos1-(yPos1[1]-yPos0[1])/2,labels=treeLabels,xpd=NA,pos=4,cex=0.7)
          
}


getAdmixFractionsTable <- function(PCset,type="PC1.1date",fitQuality,subsetPopns=NULL,popOrder=NULL,plainNames=TRUE){
  
            PCinfo = PCset[[type]]
            PCdataA <- PCinfo$a[,-1]
            PCdataA <- PCdataA[,gsub("LOO-SP_PR_HM_NA.","",rownames(fitQuality))]
            if(!is.null(subsetPopns)) PCdataA <- PCdataA[,colnames(PCdataA)%in%subsetPopns]
            if(!is.null(popOrder)) PCdataA <- PCdataA[,popOrder]
            
            #PCdataA <- PCdataA[,order(match(gsub('.null','',colnames(PCdataA)),ord))]
            xA <- apply(as.matrix(PCdataA),2,function(x) x[-length(x)]*x[length(x)])
            
            PCdataB <- PCinfo$b[,-1]
            PCdataB <- PCdataB[,colnames(PCdataA)]
            xB <- apply(as.matrix(PCdataB),2,function(x) x[-length(x)]*x[length(x)])
            
            dates <- fitQuality[match(paste("LOO-SP_PR_HM_NA.",colnames(PCdataA),sep=""),rownames(fitQuality)),]
            
            names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[3,]
            if(plainNames==T) names <- sapply(gsub('.null','',colnames(PCdataA)),getPopnLabel,popnLabels=popnLabels)[1,]
            names[grep("null",colnames(PCdataA))] <- paste(names[grep("null",colnames(PCdataA))],".null",sep="")
            
            theseA = rowSums(xA)>0 
            theseB = rowSums(xB)>0
            allDonorsA = names(which(apply(PCdataA[-nrow(PCdataA),],1,function(x) sum(x>tol)>0)))       # at least one target above tolerance
            allDonorsB = names(which(apply(PCdataB[-nrow(PCdataB),],1,function(x) sum(x>tol)>0))) 
          
            donorOrderA = allDonorsA[order(rowSums(xA[allDonorsA,],na.rm=T))]
            donorOrderB = rev(allDonorsB[order(rowSums(xB[allDonorsB,],na.rm=T))])
          
            labsA = sapply(paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderA)),getPopnLabel,popnLabels)[1,]
            labsB = sapply(paste0("donorPopn_",gsub("[^[:digit:]]","",donorOrderB)),getPopnLabel,popnLabels)[1,]
            
            dataA = round( t(PCdataA[donorOrderA,]),4)
            dataB = round( t(PCdataB[donorOrderB,]),4)
            colnames(dataA)= labsA
            colnames(dataB)= labsB
            rownames(dataA)= rownames(dataB) = names
            
            totalA = t(PCdataA["totalPCprop",])
            totalB = t(PCdataB["totalPCprop",])
            colnames(totalA) = "Total fraction (minor)"
            colnames(totalB) = "Total fraction (major)"
            toPrint = cbind(dataA,totalA,dataB,totalB)
            
        return(toPrint)           
}


getPopnCentres <- function(pops,split.matrix,stuff,fun="median",columns=c("Xmuni.ave.grand.precise", "Ymuni.ave.grand.precise")){
  
  centres <- data.frame(cbind(rep(NA,length(pops)),rep(NA,length(pops))))
  rownames(centres) <- pops
  colnames(centres) <- c("X","Y")
  
  for(pop in pops){
    print(pop)
    npop <- as.numeric(str_split(pop,"_",2)[[1]][2])
    inds <- rownames(split.matrix)[split.matrix[,2]==npop]
    print(inds)
    coords <- apply(stuff$D[inds,columns],2,function(x) do.call(fun,args=list(x,na.rm=T)))
    centres[pop,] <- coords
    
  }
  
  SpainGeogOrder=rownames(centres)[order(centres$X)]
  
  centres$map.name <- rownames(centres)
  centres[is.na(centres[,1]),c("X","Y")] <- c(0,0)
  
  return(centres)
}

#######################
# functions for matrix manipulation

sim2dist <- function(mx) {
 s1 = matrix(1,length(diag(mx)),1)%*%t(diag(mx))
 out = s1 + t(s1) - 2*mx
 rownames(out) = rownames(mx)
 colnames(out) = colnames(mx)
 return(out)
  
}


#######################
# function to rotate dendro objects
flip_dendro = function(dend, direction="x"){
    newdend=dend
    if("x"%in%direction){
      #newdend$segments[,c("x","xend")] = -dend$segments[,c("x","xend")] + 
      #                                  max(dend$segments[,c("x","xend")]) + 
      #                                  min(dend$segments[,c("x","xend")])
      newdend$segments[,c("x","xend")] = -dend$segments[,c("x","xend")] + 
                                         + dim(dend$labels)[1] + 1
                                              
      newdend$labels[,3] = rev(dend$labels[,3])
    }
    if("y"%in%direction){
      newdend$segments[,c("y","yend")] = -dend$segments[,c("y","yend")] + max(dend$segments[,c("y","yend")]) + min(dend$segments[,c("x","xend")])
    }
    
  return(newdend)
}

#######################
# new version of plotFinestructure (with rotated matrix and dendrogram)


plotFinestructure2<-function(tmpmat,# this is the heatmap to be plotted
  	labelsx,# the x,y labels to be plotted
		labelsy=NULL,# assumed equal to labelsx unless otherwise specified
		labelsatx=NULL, # location of the X labels - note that the heatmap X ranges 1:dim(tmpmat)[1] and good locations are (1:dim(tmpmat)[1]) -0.5
		labelsaty=NULL, # assumed equal to labelsatx unless otherwise specified
		cols=NULL,# colour scale. Will be generated by MakeColors if NULL
		dend=NULL, # optional dendrogram to be placed at the side
    extendDend=NULL, # extend dendro leaves by this fraction of the dendro height (only valid if class(dend)="dendro")                             
    dendSide = "left",
		optpts=NULL, # optional points to add, e.g. MAP state
		labmargin=8, # margin for x,y labels
		layoutd=0.2, # proportion of the plot to be used by the dendrogram
		layoutf=0.1,# proportion of the plot to be used by the scale
		cex.axis=0.5, # cex.axis applies to the x,y labels
		crt=0, # character rotation of x,y labels (x is rotated the opposite direction)
		colscale=NULL, # specify a colour scale (default: ???)
		text.colx=NULL, # colour of the labels (can be a vector, is recycled
    text.coly=NULL,
    text.col=NULL,
		ignorebelow=0, # if text.col=0, then we white out everything below ignorebelow
		nodePar=list(cex=0,lab.cex=0.5,las=2),# nodePar as seen by plot.dendrogram
		edgePar=list(p.col="transparent",p.border="transparent",p.lwd=1,t.srt=90,t.off=-0.1,t.cex=0.8,t.col="transparent"),# edgePar as seen by plot.dendrogram
		dendmar=c(0,1,1,0), # additional modification to the margins of the dendrogram
		scalemar=c(0,2,5,2), # additional modification to the margins of the scale
		hmmar=c(0,0,1,0), # additional modification to the margins of the heatmap
		cex.scale=1, # size of the characters in the scale
		scalelocs=NULL, # optional firced positions of the scale
		scalenum=10, # if scalelocs=NULL, the number of points to label on the scale
		scalesignif=3, # number of significant digits on scale
		scalelabel="", # label for the scale
		optcex=1.0, # size of the optional points
		optpch=20, # pch for the optional points
		optcol="white", # colour for the optional points
		labelsoff=c(1.2,1.2), # "margin" like, distance from colour/rotated labels and the axis
		tickmarks=1,  # size of tickmarks when using colour/rotated labels
    col.ticks="black",  # colour of tickmakrs
    main=NULL)
{
 # png(filename,height=2000,width=2200,res=200)
  
# rotate matrix 90 degrees and reverse dendrogram (i.e reverse the leaf order)
tmpmat=rotate90(tmpmat)

# extend dendro leaves
if(class(dend)=="dendro"){
  if(!is.null(extendDend)) {
    extendDend2 = extendDend*max(dend$segments$y)
      
    dend$segments$yend[dend$segments$yend==0] = -extendDend2
    dend$segments$yend = dend$segments$yend + extendDend2
    dend$segments$y = dend$segments$y + extendDend2
  }  
}

if(dendSide =="left"){
  if(class(dend)=="dendro") {
      dend = flip_dendro(flip_dendro(dend,"x"),"y")
  } else {
    dend=rev(dend)
  }

layout(matrix(c(1,3,2,4), 2, 2, byrow=TRUE),widths=c(layoutd,1-layoutd),height=c(layoutf,1-layoutf))
}

if(dendSide =="right"){
  if(class(dend)=="dendro") {
      dend = flip_dendro(dend,"x")
  } else {
    dend=rev(dend)
  }
  
layout(matrix(c(3,1,4,2), 2, 2, byrow=TRUE),widths=c(1-layoutd,layoutd),height=c(layoutf,1-layoutf))
  
}

#layout(matrix(c(3,1,4,2), 2, 2, byrow=TRUE),widths=c(layoutd,1-layoutd),height=c(1-layoutf,layoutf))

if(is.null(cols)) {
  #cols<-MakeColorYRP()
  #cols <- c(colorRampPalette(c("yellow","red","blue"),interpolate="linear")(100),"black") # as above, but with a black final for capped values
  cols = c(colorRampPalette(c("white","yellow","red","blue","#00005C"),interpolate="linear")(100),"black")
}
  
## 1. TOP LEFT
par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),type="n",axes=F,xlab="",ylab="")# null plot for top right

## 2. BOTTOM LEFT: DENDROGRAM
par(mar=c(labmargin,0,0,0)+dendmar+ c(0,0,hmmar[3],0))

if(is.null(dend)){
	plot(c(0,1),c(0,1),type="n",axes=F,xlab="",ylab="") # null plot for top right
}else {
  if(class(dend)=="dendro"){    
    print("plotting dendro object (rotated)")
      #print(c(min(dend$segments$y),max(dend$segments$yend)))
    
    plot(0, xlim=c(min(dend$segments[,c("y","yend")]),max(dend$segments[,c("y","yend")])),
         ylim=c(0.5,dim(dend$labels)[1]+0.5),axes=F,xlab="",ylab="",type="n",yaxs = "i")
    
    segments(x0=dend$segments$y,y0=dend$segments$x,x1=dend$segments$yend,y1=dend$segments$xend)
        
	} else {
    print("plotting dendrogram object (rotated)")
    if(dendSide=="right") TreeSide = T else TreeSide=F
    plot_horiz.dendrogram(dend,horiz=T,axes=FALSE,yaxs = "i", leaflab = "none",nodePar=nodePar,edgePar=edgePar,side=TreeSide)
	}
}
if(is.null(colscale)) colscale<-c(max(ignorebelow,min(tmpmat[tmpmat>ignorebelow])),max(tmpmat))

## 3. TOP RIGHT: SCALE
if(dendSide == "left") par(mar=c(0,0,0,labmargin)+scalemar)
if(dendSide == "right") par(mar=c(0,labmargin,0,0)+scalemar)

colindex<-matrix(seq(min(colscale),max(colscale),length.out=100),ncol=1,nrow=100) # colour scale
image(1:100,1,colindex,xaxt="n",yaxt="n",xlab=scalelabel,ylab="",col=cols,zlim=range(colindex),main=main)
if(is.null(scalelocs)){ # fill the range with scale labels
	scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=scalenum)
	scalephysicalpos <- seq(1,100,length.out=scalenum)
}else{
	stop("NOT IMPLEMENTED...")
}
axis(3,at=scalephysicalpos,labels=signif(scalelocs,scalesignif),las=1,cex.axis=cex.scale,padj=1)

## 4. BOTTOM RIGHT: MAIN HEATMAP
if(dendSide == "left") par(mar=c(labmargin,0,0,labmargin)+hmmar)
if(dendSide == "right") par(mar=c(labmargin,labmargin,0,0)+hmmar)

if(is.null(labelsaty)) { # labels y positions (columns on matrix)
	labelsaty<-seq(1:dim(tmpmat)[2])
}
if(is.null(labelsatx)) { # labels x positions (rows on matrix)
	labelsatx<-seq(1:dim(tmpmat)[1])
}
if(is.null(labelsy)) { # y labels
	labelsy<-rev(labelsx)
}

#plot the heatmap

image(1:dim(tmpmat)[1],1:dim(tmpmat)[2],tmpmat,xaxt="n",yaxt="n",xlab="",ylab="",col=cols,zlim=colscale)
# draw the axes
if(crt==0 && is.null(text.colx)) {
axis(1,at=labelsatx,labels=labelsx,las=2,cex.axis=cex.axis,col.ticks=col.ticks)
if(dendSide == "left") axis(4,at=labelsaty,labels=labelsy,las=2,cex.axis=cex.axis,col.ticks =col.ticks)
if(dendSide == "right") axis(2,at=labelsaty,labels=labelsy,las=2,cex.axis=cex.axis,col.ticks =col.ticks)

print("using defaults")
}else{
    text(labelsatx, -labelsoff[1], srt = 90-crt, adj = 1,
          labels = labelsx, xpd = TRUE,cex=cex.axis,col=text.colx)
    if(dendSide == "left") text(max(labelsatx)+labelsoff[2], labelsaty, srt = crt, pos=4, offset = 1,
          labels = labelsy, xpd = TRUE,cex=cex.axis,col=text.coly)
    if(dendSide == "right") text(min(labelsatx)+labelsoff[2], labelsaty, srt = crt, pos=2, offset = 1,
          labels = labelsy, xpd = TRUE,cex=cex.axis,col=text.coly)
    
    print("using labelsoff")
	if(tickmarks>0) {
    print('plotting tickmarks')
		axis(1,at=labelsatx,labels=rep("",length(labelsatx)),las=2,cex.axis=cex.axis,col.ticks=col.ticks)
		if(dendSide == "left") axis(4,at=labelsaty,labels=rep("",length(labelsaty)),las=2,cex.axis=cex.axis,col.ticks=col.ticks)
    if(dendSide == "right") axis(2,at=labelsaty,labels=rep("",length(labelsaty)),las=2,cex.axis=cex.axis,col.ticks=col.ticks)
	}

}
## add optional points
if(!is.null(optpts))tmp<-sapply(1:dim(optpts)[1],function(x){
  plist<-which(optpts[x,]==1);
  xpos<-rep(x,length(plist))
  ypos<-plist
#  if(yrev) ypos<-dim(tmpmat)[1] - plist + 1
  points(xpos,ypos,pch=optpch,cex=optcex,col=optcol)
  invisible(NULL)
})

## DONE

}



###############
# get set of rectangles for tree (given a level of the tree)

getRects <- function(v,thickness=10){
  # v is a vector of contiguous labels like that input into 'niceLabels'
  rects=rle(v)$values
  nRects=length(rects)
  heights=cumsum(rle(v)$lengths)
  xl=rep(-thickness,nRects)
  yb=c(0,heights[-nRects])
  xr=rep(0,nRects)
  yt=heights
  out=list(xl,yb,xr,yt)
}


###############
# plot points/rectangles on side of finestructure heat plot (haven't actually used this I don't think!)
#pointOrder = paste("SpainPopn_",SPAIN.split.matrix[labelOrder,2],sep="")
#rectOrder = paste("SpainPopn_",new.split.matrix.A2_v1_v7a[labelOrder,14],sep="")

plotFineColours <- function(pointOrder,rectOrder,boundaries=T){
  
if(is.null(pointOrder)) splitBoundaries <- cumsum(rle(rev(rectOrder))$lengths) else splitBoundaries <- cumsum(rle(rev(pointOrder))$lengths)

# shapes and colours for points
if(!is.null(pointOrder)){
    
  l=unique(pointOrder)  # l is list of populations in order of their plotting
  
  shapes=recieverLegend2[l,"shape"]
  cols=recieverLegend2[l,"colour"]
  
  
  labPosR <- splitBoundaries - rle(rev(pointOrder))$lengths/2
  labPosC <- cumsum(rle(pointOrder)$lengths) - rle(pointOrder)$lengths/2
  
  pointBorderCols = rep("black",length(labPosC))
  pointBorderCols[l=="SpainPopn_2"] = "gray"
  
}

# colours for rectangles
collapsed = rle(rectOrder)$values
labcols=recieverLegend2[rev(collapsed),"colour"]

treeRects = getRects(rev(rectOrder),thickness=40)
treeRectsBottom = getRects(rectOrder,thickness=40)
      
  

# run plotting
      if(boundaries) abline(h=splitBoundaries)
      par(xpd=NA)
      # right rectangles and points
      rect(length(rectOrder) - treeRects[[1]],treeRects[[2]],length(rectOrder) + treeRects[[3]] + 0.5,treeRects[[4]],col=labcols,border=NA) +      
      if(!is.null(pointOrder)) points(rep(length(pointOrder)+20,length(labPosR)),labPosR,bg=rev(cols),pch=rev(shapes),col=rev(pointBorderCols),xpd=NA)
      # bottom rectangles and points  
      rect(treeRectsBottom[[2]],treeRectsBottom[[1]],treeRectsBottom[[4]],treeRectsBottom[[3]] - 0.5,col=rev(labcols),border=NA) +            
      if(!is.null(pointOrder)) points(labPosC,rep(-20,length(labPosC)),bg=cols,pch=shapes,col=pointBorderCols,xpd=NA)      
  
}


plotFineColours2 <- function(pointOrder,rectOrder,boundaries=T,side="left",thickness=40,
                             nWidth=NULL,buffer=0.5,recieverLegend2=recieverLegend2){
  
if(is.null(nWidth)) nWidth=length(rectOrder) # by default the x-width is the number of samples in rectOrder. This won't be correct if the matrix is not square!

if(is.null(pointOrder)) splitBoundaries <- cumsum(rle(rev(rectOrder))$lengths) else splitBoundaries <- cumsum(rle(rev(pointOrder))$lengths)

# shapes and colours for points
if(!is.null(pointOrder)){
    
  l=unique(pointOrder)  # l is list of populations in order of their plotting  
  shapes=recieverLegend2[l,"shape"]
  cols=recieverLegend2[l,"colour"]
  
  labPosR <- splitBoundaries - rle(rev(pointOrder))$lengths/2 + 0.5
  labPosC <- cumsum(rle(pointOrder)$lengths) - rle(pointOrder)$lengths/2  + 0.5
  
  pointBorderCols = rep("black",length(labPosC))
  values = rgb2hsv(col2rgb(cols))[3,]
  pointBorderCols[values < 0.2] = "gray" # dark colours get gray border
  
}

# colours for rectangles
collapsed = rle(rectOrder)$values
labcols=recieverLegend2[collapsed,"colour"]

if(side%in%c("right","left")) {
  treeRects = getRects(rev(rectOrder),thickness=thickness)
  
}
if(side=="bottom") {
  treeRects = getRects(rectOrder,thickness=thickness)
  
}
      

# run plotting
      if(boundaries) abline(h=splitBoundaries)
      par(xpd=NA)
      
      if(side=="right"){
        # right rectangles and points
        xL = nWidth + treeRects[[3]] + buffer + 0.5
        xR = nWidth - treeRects[[1]] + buffer + 0.5
        rect(xL,treeRects[[2]] + 0.5,xR,treeRects[[4]] + 0.5,col=rev(labcols),border=NA) + 
        if(!is.null(pointOrder)) {
          pR = nWidth + treeRects[[3]][1] + buffer + thickness/2 + 0.5
          points(rep(pR,length(labPosR)),labPosR,bg=rev(cols),pch=rev(shapes),col=rev(pointBorderCols),xpd=NA)
        }
      }

      if(side=="left"){
        # left rectangles and points
        xL = treeRects[[1]] - buffer + 0.5
        xR = treeRects[[3]] - buffer + 0.5
        rect(xL,treeRects[[2]] + 0.5,xR,treeRects[[4]] + 0.5,col=rev(labcols),border=NA) +      
        if(!is.null(pointOrder)) {
          pL = treeRects[[1]][1] - buffer + thickness/2 + 0.5
          points(rep(pL,length(labPosR)),labPosR,bg=rev(cols),pch=rev(shapes),col=rev(pointBorderCols),xpd=NA)
        }
      }

      if(side=="bottom"){
        # bottom rectangles and points  
        rect(treeRectsBottom[[2]] + 0.5,treeRectsBottom[[1]],treeRectsBottom[[4]] + 0.5,treeRectsBottom[[3]] - buffer,col=labcols,border=NA) +            
        if(!is.null(pointOrder)) points(labPosC,rep(-20,length(labPosC)),bg=cols,pch=shapes,col=pointBorderCols,xpd=NA)      
      }
}

##############################
# Function to plot all the necessry bits and pieces for a nicely labelled coancestry plot

# NOTE: the dendrogram g defines the set of individuals, and the order that will be plotted in.
# This only works for square versions!
plotNiceCoancestry <- function(g,datamatrix,SPAIN.split.matrix,SPAINpopnLabels,recieverLegend2,treeDataIndex2=NULL,
                               recThickness=40,LineWidth=1,cexPoints=1,title="Coancestry (cM)",showBoundaryLines=TRUE,
                               cex.axis=0.8,labmargin=8,pdf=FALSE){
      
    treeOrder <- as.character(g$labels[,3])
    treeDataIndex <- match(treeOrder,rownames(datamatrix))
    if(is.null(treeDataIndex2)) treeDataIndex2 = treeDataIndex
    labNames=sapply(paste("SpainPopn_",SPAIN.split.matrix[treeOrder,2],sep=""),getPopnLabel,SPAINpopnLabels)[1,]
    labels <- niceLabels(labNames)
    labels = gsub("Spain_","",labels)
    
    #### empty labels
    #labels = ""
    textcol = "black"
    splitBoundaries <- cumsum(rle(SPAIN.split.matrix[rev(treeOrder),2])$lengths)
    ticksR=cumsum(rle(rev(labNames))$lengths)
    ticksC=cumsum(rle(labNames)$lengths)
    
    splitBoundariesExtraC <- cumsum(rle(labNames)$lengths)
    
    labPosR <- splitBoundaries - rle(SPAIN.split.matrix[rev(treeOrder),2])$lengths/2
    labPosC <- cumsum(rle(SPAIN.split.matrix[treeOrder,2])$lengths) - rle(SPAIN.split.matrix[treeOrder,2])$lengths/2
    
    labcols = recieverLegend2[rev(paste0("SpainPopn_",SPAIN.split.matrix[treeOrder,2])),"colour"]
    labshapes = recieverLegend2[rev(paste0("SpainPopn_",SPAIN.split.matrix[treeOrder,2])),"shape"]
    labColShapes = niceLabels(paste(labcols,labshapes,sep="___"))
    labColShapes = labColShapes[labColShapes!=""]
    labcols = sapply(labColShapes,function(x) str_split(x,"___")[[1]][1])
    labshapes = as.numeric(sapply(labColShapes,function(x) str_split(x,"___")[[1]][2]))
    pointBorderCols = rep("black",length(labPosC))
    # give white borders to dark background colours
    dark = rgb2hsv(col2rgb(labcols))["v",]<0.56
    pointBorderCols[dark] = "gray"
    
    treeCols = niceLabels(recieverLegend2[rev(paste0("SpainPopn_",SPAIN.split.matrix[treeOrder,1])),"colour"]);
    treeCols = treeCols[treeCols!=""]
    treeRects = getRects(rev(paste0("SpainPopn_",SPAIN.split.matrix[treeOrder,1])),thickness=recThickness)
    treeRectsBottom = getRects(paste0("SpainPopn_",SPAIN.split.matrix[treeOrder,1]),thickness=recThickness)
  
    if(!pdf) png(filename,height=4500,width=5000,res=450)
    #png(filename,height=2000,width=2200,res=200)
    if(pdf) pdf(gsub(".png",".pdf",filename),height=4500/450,width=5000/450)
   p <- plotFinestructure2(datamatrix[treeDataIndex,treeDataIndex],labelsx=paste(labels,"   "),dend=g,tickmarks=1,text.colx=textcol,text.coly=textcol,labelsoff=c(10,10),
                             cols=some.colorsEnd2,cex.axis=cex.axis,cex.scale=0.8,edgePar=list(p.col="transparent",p.border="transparent",p.lwd=1,t.srt=90,t.off=-0.1,t.cex=0.8,t.col="transparent",lwd = 1.2),
                              col.ticks="transparent",dendSide="right",dendmar = c(0,0.5,0,1),labmargin=labmargin,extendDend=0.05,
                             #main=paste(prefix,": ",matrixName,"; chromV",chromvSpain,sep="")) + 
                              main=title)
   if(showBoundaryLines) p = p + abline(h=splitBoundaries,lwd=LineWidth)
      
      # extend lines beyond plot to demarkate clust + add labels and colours to dendrogram   
      par(xpd=NA)
      #plot(1:length(treeDataIndex),1:length(treeDataIndex))
      for(l in splitBoundariesExtraC) lines(y=c(-length(treeDataIndex)/50,0),x=rep(l,2),lwd=LineWidth) # extend beyond plot boundaries
      for(l in splitBoundariesExtraC) lines(x=c(-length(treeDataIndex)/50,0),y=rep(length(treeDataIndex)-l,2),lwd=LineWidth) # extend beyond plot boundaries
      
      # right rectangles and points
      pointBuffer = length(treeDataIndex)/70
      
      rect(length(treeDataIndex) - treeRects[[1]],treeRects[[2]],length(treeDataIndex) + treeRects[[3]] + 0.5,treeRects[[4]],col=treeCols,border=NA) +      
      points(rep(length(treeDataIndex)+pointBuffer,length(labPosR)),labPosR,bg=labcols,pch=labshapes,col=pointBorderCols,xpd=NA,cex=cexPoints)
      
      # bottom rectangles and points  
      #rect(treeRectsBottom[[2]],treeRectsBottom[[1]],treeRectsBottom[[4]],treeRectsBottom[[3]] - 0.5,col=rev(labcols),border=NA) +            
      points(labPosC,rep(-pointBuffer,length(labPosC)),bg=rev(labcols),pch=rev(labshapes),col="black",xpd=NA,cex=cexPoints)
      
      # left points
      points(rep(-pointBuffer,length(labPosR)),labPosR,bg=labcols,pch=labshapes,col="black",xpd=NA,cex=cexPoints)          
      par(opar)
  
      dev.off()
      
      
    return(NULL)
}


#######################
# Test of cluster assignments in two finestructure runs (e.g in SPAIN.A2.Galicia.treeplots.R)
# How often does each pair of samples end up in the same cluster (or different cluster)?

fracDiff = function(X){      
      a = t(X)[lower.tri(t(X))]  # upper triangle of X
      b = X[lower.tri(X)]  # lower triangle of X
      raw = sum(a!=b)/((dim(X)[1]^2-dim(X)[1])/2)
      rawA = sum(a[a==1]!=b[a==1])/sum(a==1)  # pairs co-clustered in A only
      rawB = sum(a[b==1]!=b[b==1])/sum(b==1)  # pairs co-clustered in B only
      return(c(raw,rawA,rawB))
    }

isDiff = function(X){
      a = t(X)[lower.tri(t(X))]  # upper triangle of X
      b = X[lower.tri(X)]  # lower triangle of X
      out = X
      out[lower.tri(out)] = a==b
      diag(out) = rep(1,dim(out)[1])
      out[upper.tri(out)] = t(out)[upper.tri(t(out))]
      
      return(out)
}

getBoth <- function(matA,matB,nClusters){
  # take two cluster assignment vectors and return a 0,1 indicator matrix showing pairwise co-assignments for each cluster vector.
  
    columnA = which(apply(matA,2,function(x) length(unique(x)))==nClusters)
    columnB = which(apply(matB,2,function(x) length(unique(x)))==nClusters)
    
    classA = matA[,columnA[1]]
    classB = matB[,columnB[1]]
    
    assignment = function(X){
      out = sapply(1:length(X), function(i) {
            a = sapply(1:length(X),function(j){
                sum(X[i]==X[j])
            })
        })
    output <- matrix(unlist(out), ncol = length(X), byrow = T)
    }
    
    clustAssignA = assignment(classA)
    rownames(clustAssignA) = colnames(clustAssignA) = rownames(matA)
        
    clustAssignB = assignment(classB)
    rownames(clustAssignB) = colnames(clustAssignB) = rownames(matB)
    
    #clustAssignB = clustAssignB[c(1,4,10,70,90,102,155,221),c(1,4,10,70,90,102,155,221)]
    #clustAssignA = clustAssignA[c(1,4,10,70,90,102,155,221),c(1,4,10,70,90,102,155,221)]
    #image(clustAssignB)
    #image(clustAssignA)
    
    # order by A
    ord = rownames(mata)[rownames(mata)%in%rownames(clustAssignA)]
    
    both=clustAssignA[ord,ord]
      B =  clustAssignB[ord,ord]
    both[lower.tri(both)] = B[lower.tri(B)]
        
    return(both)
}
    
getFracDiff <- function(versionA,versionB,nClusters=10,maxCluster=10,plotImage=F,byInd=F){
    
  print(paste(versionA,versionB,nClusters,sep="; "))
  
    mata = get(paste('new.split.matrix.',pre,'_',versionA,sep=""))
    matA = mata[!rownames(mata)%in%c(samplesToRemove,"nonGalicia"),]
    
    # remove poor samples
    matb = get(paste('new.split.matrix.',pre,'_',versionB,sep=""))
    matB = matb[!rownames(matb)%in%c(samplesToRemove,"nonGalicia"),]
    
    # subset by smallest set (or reorder matB so that it matches matA)
    if(dim(matA)[1]>dim(matB)[1]) matA = matA[rownames(matB),] else matB = matB[rownames(matA),]
        
    if(byInd){
      print('computing co-assignment matrices for each cluster')
      # for each cluster level, get matrix of pairs indicating whether the co-assignment is the same in both runs                
      bothMats = sapply(2:maxCluster,function(i) {
                both=getBoth(matA,matB,i)
                fracMat = isDiff(both)
      },simplify=F)        
      
      if(plotImage){
        image(bothMats[[3]])            
      }
      return(bothMats)
    
    } else {
            
        both=getBoth(matA,matB,nClusters)
        
        # fraction of pairs with the same co-assignment in both runs
        frac = fracDiff(both)
          
      if(plotImage){
        image(clustAssignA)    
        image(clustAssignB)    
        image(both)    
      }
          
      return(list(frac,both))
    }
}

#######################
# plotting functions for globetrotter results

plotGToutput = function(type,plots=c("bars","residuals","contrasts"),
                        extra="",bootstraps=T,bootsCol="date1.est.boot",fixYaxis=90,dateCol=dateCol){

  PCinfo = PCset[[type]]
  pops = colnames(PCinfo$a)[-1]
  typeName = types[type]
  conclusions=as.character(fitQuality$conclusion)
  names(conclusions)=gsub('LOO-SP_PR_HM_NA.','',rownames(fitQuality))
    
  if("contrasts"%in%plots){
    print(paste('getting contrasts for ',type))

    contrasts = list()
    for (pop in pops){
      #print(pop)
      contrasts[[pop]] = getContrast(targetPopn=pop,copyVectors=SUMMARYmat,PCinfo,contrastType="1")
    }
    
    contrastMat = matrix(NA,nrow=length(donorOrder),ncol=length(popOrder))
    colnames(contrastMat) = popOrder
    rownames(contrastMat) = donorOrder
      
    for (pop in popOrder){  
      contrastMat[,pop] = contrasts[[pop]][rownames(contrastMat),1]
    }
    
    print(paste('plotting contrasts...'))
    
    filename=paste(plotdir,'/GTcontrasts-',version,'-',type,extra,'-%02d.png',sep="")
      toPlot = contrastMat
      plotMixtureHeat(toPlot,filename,title=paste("Mixture contrasts : ",type,sep=""),cap=F,
                      colourSet=c("red","white","blue"),Cex=1.5,
                      ylab="target popns",xlab="donors",negativeScale=T,nCols=100,width=1500,height=1500,res=100)

  }
  
  if("bars"%in%plots){
    print(paste('plotting barplots for',type))
    
    filename=paste(plotdir,'/GTbarplots-',version,'-',type,extra,'-%02d.png',sep="")
    png(filename,width=3000,height=2000,res=200)
      par(opar)
        
      PC.barplot.dates(PCinfo=PCinfo,subsetPopns=popOrder,popOrder=popOrder,lines=Lines,fullSurr=F,sources=bestMatchSources,title=paste(typeName,", best sources",sep=""),bootsCol=bootsCol,
                         dateCol=dateCol,sampSizes=sizes,fixYaxis=fixYaxis,bootstraps=bootstraps,conclusions=conclusions,mar=c(12,5,3,15),cex.axis=0.5,lwd=0.5,cex=1)
      
      PC.barplot.dates(PCinfo=PCinfo,subsetPopns=popOrder,popOrder=popOrder,lines=Lines,title=paste(typeName,", source mixture",sep=""),bootsCol=bootsCol,
                         dateCol=dateCol,sampSizes=sizes,fixYaxis=fixYaxis,bootstraps=bootstraps,conclusions=conclusions,mar=c(12,5,3,15),cex.axis=0.5,lwd=0.5,cex=1)
      
    dev.off()
    
  }
  
  if("residuals"%in%plots){
    print(paste('getting residuals for',type))
    
    doPlot=T
    residual <- PCinfo$a[-nrow(PCinfo$a),-1]
      
    for (pop in colnames(residual)){
    
      #print(pop)
          
        a <- PCinfo$a[-nrow(PCinfo$a),pop]*PCinfo$a[nrow(PCinfo$a),pop]
        b <- PCinfo$b[-nrow(PCinfo$b),pop]*PCinfo$b[nrow(PCinfo$b),pop]
        ab = a+b
        names(ab) <- rownames(PCinfo$a)[-nrow(PCinfo$a)]
        # remove self-copying (where applicable)
        
        ab <- ab[names(ab)!=gsub(".null","",pop)]
        
        copyVectors=SUMMARYmat
        rownames(copyVectors) = sapply(rownames(copyVectors),getPopnLabel,popnLabels)[2,]
        colnames(copyVectors) = sapply(colnames(copyVectors),getPopnLabel,popnLabels)[2,]
        
        # make portugal index in copyVectors match that in PCinfo
        rownames(copyVectors)[grep("Portugal",rownames(copyVectors))] = gsub(".null","",colnames(PCinfo$a)[grep("Portugal",colnames(PCinfo$a))][1])
  
        copyVect <- copyVectors[,colnames(copyVectors)!=gsub(".null","",pop)]
        copyVect2 <- copyVect/rowSums(copyVect)
        copyVect3 = copyVect2
        rownames(copyVect3) = sapply(rownames(copyVect3),getPopnLabel,popnLabels)[2,]
      
        admixVect <- t(copyVect3[names(ab),])%*%ab
      
      if((pop%in%c("Portugal_6","Portugal_6.null"))&(version%in%c("v6","v7"))) pop2="Portugal_49" else pop2=getPopnLabel(gsub(".null","",pop),popnLabels)[2]
        
        resid = admixVect - copyVect2[gsub(".null","",pop2),]
        colnames(resid) <- pop
        rownames(resid) = sapply(rownames(resid),getPopnLabel,popnLabels)[2,]
                
        residual[,pop] <- NA
        residual[rownames(resid),pop] <- resid
    
      }
    
    forPlot <- residual#[,!is.na(residual[1,])]
    #forPlot <- forPlot[,ord]  
    filename=paste(plotdir,'/GTresiduals-',version,'-',type,extra,'.png',sep="")
    
    if(doPlot==T){
      print(paste('plotting residuals...'))
    
    plotMixtureHeat(toPlot=forPlot[,popOrder],filename=filename,
                    title="Admixture inference residuals: inferred - actual",cap=F,colourSet=c("red","white","blue"),
                    Cex=1.5,xlab="donor popn",ylab="target popn",width=1500,height=1500,res=100)
    }
    TVD = 0.5*colSums(abs(residual))
    assign(paste(type,".TVD",sep=""),TVD)
  }
  
}

plotCurveFits = function(curveFitDataAll,dateType="gen.fit.1date",fitType="intercept.fit"){

 for(pop in names(curveFitDataAll)){
  
  print(pop)
  matrix=curveFitDataAll[[pop]][[fitType]][[dateType]]
  nonNA=names(which(rowSums(!is.na(matrix))!=0))
  
  if(grepl("null",pop)) popName=paste(getPopnLabel(gsub(".null","",pop),popnLabels)[2],".null",sep="") else popName=getPopnLabel(pop,popnLabels)[2]
  
  toPlot=matrix[nonNA,nonNA]
  toPlot = toPlot[order(rowSums(toPlot)),order(rowSums(toPlot))]
  File=paste(plotdir,'/curveFits/GTcurveFits-',version,'-',dateType,'-',fitType,'-',pop,'-%02d.png',sep="")
  
  if(fitType=="rsquared.date.fit") {
     plotMixtureHeat(toPlot=toPlot,filename=File,negativeScale=F,fixedRange=c(0,1),
                    title=paste(popName,":  R^2 fit for ",dateType,sep=""),cap=F,
                    Cex=1.5,xlab="surrogate popn",ylab="surrogate popn",width=1500,height=1500,res=100)

  } else {
    plotMixtureHeat(toPlot=toPlot,filename=File,
                    title=paste(popName,": curve intercepts for ",dateType,sep=""),cap=F,
                    Cex=1.5,xlab="surrogate popn",ylab="surrogate popn",width=1500,height=1500,res=100)
  }
  
    }
}

plotCurves = function(curveFitDataAll,targetPopn,filename,dates="1date",separatePages=FALSE,plot.cex=1,...){
  
    blue.col=col2rgb("cornflowerblue")
     red.col=col2rgb(2)
     green.col=col2rgb(3)
  
    dat = curveFitDataAll[[targetPopn]]$curves
        
    means = dat[dat$V3=="scaled.data",4:ncol(dat)]
    meansNames = dat[dat$V3=="scaled.data",1:2]
    times = dat[1,4:ncol(dat)]
    fitted = dat[dat$V3=="source.fit.1date",4:ncol(dat)]
    date1 =  dat[dat$V3=="gen.fit.1date",4:ncol(dat)]
    dates2 = dat[dat$V3=="gen.fit.2date",4:ncol(dat)]
    
    plot.subtitle = getPopnLabel(gsub(".null","",targetPopn),popnLabels)[1]
    print(plot.subtitle)
    
    refit=F
  pdf(filename,bg="transparent",...)
    
  if(separatePages) par(mfrow=c(1,1),cex=plot.cex) else par(mfrow=c(2,2),cex=plot.cex) 
  
    for(i in 1:nrow(means)){
      
      pop.comb.i=meansNames[i,]
                            
      ## GENERATE FITTED CURVES:
      if(refit){
        temp=lm(means[i,]~pred)
      	temp2=lm(means[i,]~pred2)
      	news=temp$coeff
      	news2=temp2$coeff
      	predline=news[1]+pred*news[2]
      	predline.fitted=news[1]+pred*c(t(intercept.mat.fitted))[i]
      	predline.2date=news2[1]+pred2[,1]*news2[2]+pred2[,2]*news2[3]
        
      } else {
        predline.fitted = fitted[i,]
        predline.2date = dates2[i,]
        predline = date1[i,]
      }        
      xlim.plot=c(min(times),max(times))
      
      
             
      plotTitle = paste(getPopnLabel(paste0("donorPopn_",gsub("[^[:digit:]]","",pop.comb.i[1])),popnLabels)[1],
                        getPopnLabel(paste0("donorPopn_",gsub("[^[:digit:]]","",pop.comb.i[2])),popnLabels)[1],
                                     sep='--')
      
      plot(t(times),t(means[i,]),type="l",lwd=3,main=plotTitle,sub=plot.subtitle,xlab="distance (cM)",ylab="Relative probability",xlim=xlim.plot)
      	#lines(times,predline.fitted,col=rgb(blue.col[1],blue.col[2],blue.col[3],maxColorValue=255,alpha=155),lwd=4)
      	if(dates=="2dates") lines(times,predline.2date,col=rgb(red.col[1],red.col[2],red.col[3],maxColorValue=255,alpha=155),lwd=4)
      	lines(times,predline,col=rgb(green.col[1],green.col[2],green.col[3],maxColorValue=255,alpha=155),lwd=4)
      	abline(h=1,col=1,lty="dotted")
      #if (i==1) legend(max(xlim.plot),max(means[i,]),legend=c("data","1-date fitted","1-date (source)","2-date"),lty=rep(1,4),lwd=rep(2,4),col=c(1,3,"cornflowerblue",2),xjust=1,yjust=1,bty='n',cex=1) 
       #legend(max(xlim.plot),max(means[i,]),legend=c("data","1-date fitted","2-dates fitted"),lty=rep(1,4),lwd=rep(2,4),col=c(1,3,2),xjust=1,yjust=1,bty='n',cex=1) 
      if((means[i,1] - means[i,ncol(means)])>0)	legendPos="topright" else legendPos="bottomright"
      legend(legendPos,legend=c("data","1-date fitted","2-dates fitted"),lty=rep(1,4),lwd=rep(2,4),col=c(1,3,2),xjust=1,yjust=1,bty='n',cex=1) 
      
        }
   
  dev.off()
    
}

# e.g 
# plotCurves(curveFitDataAll,"SpainPopn_1",width=10,height=7)


