########################################
### Create density maps (based on continous values)  ####
#######################################
 
DIR="~/Documents/ClareDPhil/DPhil/Spain/densityMaps";#setwd(DIR)
imageHeight=1500
#source("../Other/ClaresFunctions.R")

################## Variables to be set when running this.
#dir <- paste(prefix,"_v",chromv,"_v",finev,"-",filter,"Km-clusters-PropTotalGaussianDensity_DIR",sep="")
#dir = 'testing_DIR'
#inputVector = maxCert
#colorSet = c("blue","purple","red")
#pointsPositions = D
#fixLims=c(0.5,1)
#filename="test.png"
#resolution="low" #Do you want high resolution or low resolution?
#model="Gaussian"
#reCompute=TRUE
##################

################## Build spatial grids if needed otherwise just load from file in next step ####################
#source('../densityMaps/getSpatialGrid_v2.R')
##################

File = paste0(filename,model,"-sd",SD,"-",resolution,"res.png")
if(setSDbyData) File = paste0(filename,model,"-VariableSD",SD,"-",resolution,"res.png")

# can we just load the data?
plotsFile = gsub(".png",".RData",File)

################## Get grid information and sample positions
filter = 80
if( !"distances" %in% ls() ) {
  #if(resolution=="low") load(paste0(DIR,'/Port_A2-SpatialDistances-',filter,'Km-lowRes.Rdata'),verbose=T)
  if(resolution=="low") {
    load(paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-lowRes.Rdata'),verbose=T)
    SP_grd=sp_grd2
  }
  
  if(resolution=="high"){
    load(paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-highRes.Rdata'),verbose=T)  # this is the new projected version
    distances = distancesZoom
    SP_grd=sp_grdZoom2
  } 
}

##################
#rm("distancesZoom") # remove this large object for workspace unless really needed!

################## Set Gaussian Parameters 
scaleFactor = 10000
#SD=3.5   # ----> this what was used in the final version for SPAIN A2!!!!!!!!
#SD=4
#SD=1/3.5
#SD=0.2
##################
##################
# samplesFor Gaussian
samplesForGaussian = names(inputVector)
if(!PlotPortugalToo) samplesForGaussian = samplesForGaussian[!grepl("POPRES",samplesForGaussian)]  
##################

################## subset data
#distancesSubset = distances[,sampleIDs%in%samplesForGaussian] # only include samples in the input vector.
distancesSubset = t(distances[sampleIDs%in%samplesForGaussian,]) # only include samples in the input vector.
sampleIDs2 = sampleIDs[sampleIDs%in%samplesForGaussian]
##################

################## Compute distances (i.e contributions from each individual) scaled
if( ( !basename(plotsFile)%in%list.files(dirname(plotsFile)) ) | (reCompute) ) {
  
  print("Computing distances...")

if(setSDbyData){

  # Load results from getVariableBandwidth.R. Note: the matrix d must be in the order of samples in the vector sampleIDs2
  load(paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-weights-variableBWexpBalloon-',SD,'-',resolution,'Res.Rdata'),verbose=TRUE)

#  scaledDistances = distancesSubset/scaleFactor
  ### It turns out this doesn't work very well!
  #forSD = rowSums(scaledDistances^2); forSD=forSD/min(forSD) #scale distances by smallest ones
  #forSD = rowSums(scaledDistances); forSD=forSD/min(forSD) #scale distances by smallest ones
  #newSDs = SD*forSD
  
  # Radius within which we have at least 10 points.
  # Pick bandwidth such that sum of all points (weighted by a gaussian without variance-scaled) sums to some fixed value.
  
  # Alternatively: 1. Balloon estimator: scale by the (relative) denisty of data around that point (using the original fixed bandwidth)!
  # Results file called "*variableBW"
#  X <- D[sampleIDs2,"Xmuni.ave.grand.precise"]
#  Y <- D[sampleIDs2,"Ymuni.ave.grand.precise"]

#  grand.points <- SpatialPoints(cbind(X,Y))
#  proj4string(grand.points) <- CRS(proj4string(espMapTotal)) # <-- coordinates of indiivuals come from original coordinates in espMapTotal
#  sp.grand.points <- spTransform(grand.points,CRS(proj4string(espMapTotalProj)))  # in grid coordinates
#  sp_grd2Points = sp_grd2; gridded(sp_grd2Points) = FALSE  # convert grid to points (centres of grids)
  
  #samplePoints = t(sp.grand.points@coords)/scaleFactor
  #gridPoints = sp_grd2Points@coords/scaleFactor
  #cov = sqrt(SD)*diag(2)
  
  #distances2D = sapply(1:nrow(gridPoints),function(x){t(samplePoints - gridPoints[x,])},simplify = FALSE)
  
 # dens2 = sapply(distances2D,function(x){
  #  dmvnorm(x,mean=c(0,0),sigma = cov,log=TRUE)
  #})
  
  #581769,-490775 
  
  #dens <- sapply(1:nrow(scaledDistances),FUN=function(x){
  #      dnorm(scaledDistances[x,],mean=0,sd=SD,log = TRUE)
  #  })
  
  #pointDens = 1e-10 + (colSums(exp(dens))/dim(dens)[1]) # 1e-10 is to avoid dividing by zero where density is so low
  #newSDsLog = 0.5 * ( ( (1/length(pointDens))*sum(log(pointDens)) ) - log(pointDens) )
  #newSDs = SD*exp(newSDsLog)
  
  #newSDs = SD*max(pointDens)/pointDens # SD is the minimum bandwidth, and it increases with lower sampling densities
    
  # select a sensible pilot bandwidth!
  #hpi = Hpi.diag(x=t(samplePoints))
  #cov = diag(2)*mean(diag(hpi))
  
  #KDE = kde(x = t(samplePoints),eval.points = gridPoints,H = cov)
  # new = (G/f(x))^0.5. this is the log of that.
  #newSDsLog = 0.5 * ( ( (1/length(KDE$estimate))*sum(log(KDE$estimate)) ) - log(KDE$estimate) )
  #newSDs = cov[1]*exp(newSDsLog)
  
  # apply adaptive kernel to each gridpoint
  #I = diag(2)
  #d = sapply(1:length(distances2D),function(x){
  #  dmvnorm(distances2D[[x]],mean=c(0,0),sigma = newSDs[x]*I,log=TRUE)
  #})
  #d=exp(d)
  
  # Alternatively: 2. Sample point estimator: scale the bandwidth for each sample point, rather than the estimated region.
  # Didn't seem to work...
  #pointDens = 1e-10 + (rowSums(dens)/dim(dens)[2])
  #newSDs = max(pointDens)/pointDens # SD is the minimum bandwidth, and it increases with lower sampling densities
  #d <- t( sapply(1:ncol(scaledDistances),FUN=function(x){
  #    dnorm(scaledDistances[,x],mean=0,sd=newSDs[x])
  #  }) )
  
  # Plot how the bandwidth changes!
  #sdPlot = plotDenistyContinuous(newSDs,colourSet=mapsColourSet,title="",sp_grd=sp_grd2,topAlpha=0.5,fixedLims=c(0,20),
  #transparencyVector=transparencyVector,notinthesea = Notinthesea,mainMap="SpainMap_under2Proj")
  #sdPlot[[2]]
  
  
  if(model=="t-dist") {
    d <- sapply(1:nrow(scaledDistances),FUN=function(x){
      dt(scaledDistances[x,],df=1/newSDs[x])
    })
  }
  if(model=="Gaussian") {
   # d <- sapply(1:nrow(scaledDistances),FUN=function(x){
  #    dnorm(scaledDistances[x,],mean=0,sd=newSDs[x])
  #  })
  }

} else {
  
if(model=="Gaussian") d <- apply(distancesSubset/scaleFactor,MARGIN=1,FUN=dnorm,0,SD) # not zoomed version
if(model=="t-dist") d <- apply(distancesSubset/scaleFactor,MARGIN=1,FUN=dt,SD) # not zoomed version. SD is actually degrees of freedom for the t-distribution.

}

# matrix of size n x grid-size
sumd <- colSums(d) # sum distance-based contributions across each grid-point
################## 



################## Set the output directory
if( "dir" %in% ls() ) system(paste('mkdir ',DIR,"/",dir,sep=""))
################## 

################## Calculate weighted grid colors. 
# I.e multiply all distance weights by the values in the input vector; then divide by the total distance-based contribution to that point.
if( length(inputVector[sampleIDs2])!=length(inputVector) ) print("WARNING: not all of your samples are going to be contributing to the density. Will not plot these at all.")
out = inputVector[sampleIDs2]*d
toPlot <- rowSums(t(out)/sumd)
#toPlot <- colSums(t(out)/rowSums(d))

if(addTransparency) transparencyVector = sumd else transparencyVector=NULL

if(PlotPortugalToo){
  
    spainPoints = pointsPositions[rownames(pointsPositions)%in%names(inputVector),c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")]
    portPoints = pointsPositions[rownames(pointsPositions)%in%names(inputVector),c("Xmuni.mixture","Ymuni.mixture")]
    colnames(portPoints) = colnames(spainPoints)
    myPoints = rbind(spainPoints[!is.na(spainPoints[,1]),],portPoints[grepl("POPRES",rownames(portPoints)),])
    inPort = over(makeSpPoints(myPoints),PortugalMap);
    inPort = rownames(inPort)[!is.na(inPort[,1])]
    inPort = inPort[inPort%in%names(inputVector)]
      
    # just plot the mean value for portugal-based points!
    portDensity = rep(mean(inputVector[inPort]),dim(d)[2]); # a vector of length number of grid-points
}

} else {
  load(plotsFile,verbose=TRUE)
}

##################

################## fix some parameters
if(is.null(fixedLims)) fixedLims = c(min(inputVector[sampleIDs2]), max(inputVector[sampleIDs2]))
if(is.null(title)) title =""
if(is.null(scaleLab)) scaleLab = "cM"

################## plot this!

if( ( !basename(plotsFile)%in%list.files(dirname(plotsFile)) ) | (reCompute) ){
  
  print("making density plots...")
  
  Notinthesea = which(!is.na(over(SP_grd,espMapTotalSansCanProj)[,1]))
  Plots <- plotDenistyContinuous(toPlot,colourSet=colorSet,title="",sp_grd=SP_grd,topAlpha=0.5,fixedLims=fixedLims,
                                 transparencyVector=transparencyVector,notinthesea=Notinthesea,mainMap = "SpainMap_under2Proj")
  Plots = Plots[[2]]
  
  if(PlotPortugalToo){
    p = plotDenistyContinuous(portDensity,colourSet=colorSet,title="",sp_grd=SP_grd,topAlpha=0.5,fixedLims=fixedLims,
                                 transparencyVector=transparencyVector,notinthesea=InPortOnly,mainMap = "SpainMap_under2Proj")
  
    Plots = Plots + p[[2]]
  }

    
    # save the results        
    print("density plots saved in...")
    print(plotsFile)
    save(Plots,file=plotsFile)

}


################## plot the points!
#G = update(SpainMap_under2,col.regions="white") +  Plots[[2]] + update(commapborders2,lwd=0.5,col=hsv(0,0,0.6))
#G
G = update(SpainMap_under2Proj,col.regions="white") + Plots + update(commapborders2Proj,lwd=0.5,col=hsv(0,0,0.6))

png(File,width=imageHeight*spainAspect,height=imageHeight,res=150)


plotPointsOnMaps2(D=pointsPositions[rownames(pointsPositions)%in%names(inputVector),],dend=NULL,zoom=F,split.matrix=NULL,m=2,
                              showLines=F,plotDendro=F,dendroPlot=NULL,baseMap=G,
                              levelsVector=inputVector,shapesVector=symbols,levelSplits=10,
                              extraTitle="",pointSize=1,LWD=2,minSplitSize=0,
                              title=title,
                              trans=1,darken=1,showScale=T,colourSet=colorSet,
                              fixedLims=fixedLims,histXlims=fixedLims,extraPoints=NULL,
                              mapLabels=T,newPage=T,numToColours=NULL,pointLWD=0.6,
                              histBgCol="transparent",jitter="max",morePoints=NULL,scaleLab=scaleLab,
                              histBinWidth=diff(fixedLims)/50,bbox=rbind(c(-9.5,-5),c(40,44)),pointsPCH=21)
dev.off()

