########################################
###  Function to build a spatial grid (and save it).
###  THIS VERSION FIRST PROJECTS THE MAP, THEN COMPUTES THE GRID, so each cell is equal area in real life.
#######################################

DIR="~/Documents/ClareDPhil/DPhil/Spain/densityMaps";setwd(DIR)
source("../Other/ClaresFunctions.R")

######### set filter parameters for samples selected in grid
filter = 80
#########

#########
# Project map into Distance coordinates ==> already done in ClaresFunctions. Look for myProjection object
#########
#### myProjection====> Based on centre of Madrid.
print(myProjection)

######### Get grid
# low resolution (include canarias here)
bb <- bbox(espMapTotalProj)
cs <- c(1,1)*8000 # cell size (in units of proj4string(espMapTotalProj). i.e meters)
#cs <- c(1,1)*100000 # cell size (in units of proj4string(espMapTotalProj). i.e meters)

# 1 ft = 3.28084 m;  1degree X-direction = 93606.42 meters
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd
# cellcentre.offset 923018 129964
# cellsize           19685  19685
# cells.dim              8      8

sp_grd2 <- SpatialGridDataFrame(grd,
                               data=data.frame(id=1:prod(cd)),
                               proj4string=CRS(proj4string(espMapTotalProj)))

sp_grd2Points = sp_grd2; gridded(sp_grd2Points) = FALSE  # convert grid to points (centres of grids)

# 25273 cells! (c.f 23460 with original coordinate system)

# checking that points of grids are actually the centre of the grid
#h = pointDistance(sp_grd2Points,sp_grd2Points,lonlat=FALSE,allpairs = TRUE)
#spplot(makeBboxPolygon(sp_grd2,box=bbox(h[[1]])),col.regions="transparent") + spplot(h[[1]]) + spplot(sp_grd2Points,col="black",pch=10)
#sp_grd2@data$dist = h[,1]
#spplot(sp_grd2,zcol="dist") + spplot(sp_grd2Points,col="black",pch=10)


######## High-res for zoom-in versions
bb <- bbox(espMapTotalProj[rownames(espMapTotalProj@data)!="Canarias",]) # exclude canarias to save space! but include all of portugal
bb[1,1] <-  bb[1,1] - bb[1,1]/60  # expand left side to fit all of portugal in!
#spplot(makeBboxPolygon(espMapTotalProj,bb)) + spplot(PortugalMapProj) + spplot(espMapTotalProj,zcol="dummy")
            cs <- c(1,1)*3000 # cell size (in units of proj4string(espMapTotalProj). i.e meters)
            # 1 ft = 3.28084 m
            cc <- bb[, 1] + (cs/2)  # cell offset
            cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
            grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
            grd
            #cellcentre.offset -9.25580410 41.85105398
            #cellsize           0.08333333  0.08333333
            #cells.dim         31.00000000 24.00000000
            
            sp_grdZoom2 <- SpatialGridDataFrame(grd,
                                           data=data.frame(id=1:prod(cd)),
                                           proj4string=CRS(proj4string(espMapTotalProj)))
# 174,150 cells for 2.5km! 128,154 cells for 3km. For 3km the matrix is 700712 KB = 700MB
sp_grdZoom2Points = sp_grdZoom2; gridded(sp_grdZoom2Points) = FALSE  # convert grid to points (centres of grids)
      
# select samples
sampleStuff = getFilteredPointsSpain(filePrefix="",excludePont=T,split.matrix=split.matrix,addSPAINA2=TRUE)

sampleIDs <- rownames(sampleStuff$D)
X <- sampleStuff$D[,"Xmuni.ave.grand.precise"]
Y <- sampleStuff$D[,"Ymuni.ave.grand.precise"]

grand.points <- SpatialPoints(cbind(X,Y))
proj4string(grand.points) <- CRS(proj4string(espMapTotal)) # <-- coordinates of indiivuals come from original coordinates in espMapTotal
sp.grand.points <- spTransform(grand.points,CRS(proj4string(espMapTotalProj)))  # in grid coordinates

grand.points.display <- sampleStuff$D
mypoints = makeSpPoints(grand.points.display[,c("Xmuni.mixture","Ymuni.mixture")],espMapTotal)
mappoints <- list("sp.points",mypoints,col=rep("black",nrow(grand.points.display)),pch=16,cex=0.5)

spplot(espMapTotal,zcol="Total",col.regions="gray",sp.layout=list(mappoints))

sp_grd <- spTransform(sp_grd2Points,CRS(proj4string(espMapTotal))) # grid back to long/lat coordinates for plotting ==> only if needed!
distances = pointDistance(sp.grand.points,sp_grd2Points,lonlat=FALSE)

sp_grdZoom <- spTransform(sp_grdZoom2Points,CRS(proj4string(espMapTotal)))
distancesZoom <- pointDistance(sp.grand.points,sp_grdZoom2Points,lonlat=FALSE) # this is a long step - N*Ncells distances to compute, in metres

sampleDistances <- pointDistance(sp.grand.points,sp.grand.points,lonlat=FALSE,allpairs = TRUE)

# get part of grid that's in portugal only (for high resolution only)

#notintheseaSpain <- which(!is.na(over(sp_grdZoom2,SpainWithoutPortMapProj)[,1])) 
notintheseaSpain <- which(!is.na(over(sp_grdZoom2,espMapTotalSansCanProj)[,1])) 
notintheRealSea <- which(!is.na(over(sp_grdZoom2,SpainPortugalMapProj)[,1])) 
portBox = makeBboxPolygon(PortugalMapProj,bbox(PortugalMapProj@polygons[[1]]@Polygons[[1]])+cbind(c(0,-20000),c(-40000,5000)))
# rotate: 
rot=0.2; rotated = t(t(portBox@polygons[[1]]@Polygons[[1]]@coords)-rowMeans(bbox(portBox)))%*%cbind(c(cos(rot),sin(rot)),c(-sin(rot),cos(rot)))
portBox = spPolygons(t(rowMeans(bbox(portBox))+t(rotated)))
proj4string(portBox) = proj4string(PortugalMapProj)
#spplot(portBox) + spplot(PortugalMapProj,zcol="dummy")
inPortBbox <- which(!is.na(over(sp_grdZoom2,portBox))) 
InPortOnly <- intersect(inPortBbox,notintheRealSea); InPortOnly <- InPortOnly[!InPortOnly%in%notintheseaSpain]

# these are large! (1 GB!!)
#save(distances,distancesZoom,sp_grd,sp_grd2,sp_grdZoom,sampleDistances,file=paste(prefix,'-SpatialDistances-',filter,'Km.Rdata',sep=""))
#save(sampleIDs,file=paste(prefix,'-SpatialDistances-',filter,'Km-sampleIDs.Rdata',sep=""))

save(sampleIDs,distances,sp_grd2,sampleDistances,file=paste(prefix,'-SpatialDistances-',filter,'Km-lowRes.Rdata',sep=""))
save(sampleIDs,InPortOnly,distancesZoom,sp_grdZoom2,sampleDistances,file=paste(prefix,'-SpatialDistances-',filter,'Km-highRes.Rdata',sep=""))

#################################################################################
