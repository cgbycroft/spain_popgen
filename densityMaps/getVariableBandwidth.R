######################
# Kernel smoothing
######################

# Input is a set of 2D sampling points (may have a value) + set of regular 2D grid-points.
# Output is a set of weights (d) for each gridpoint and sample point pair, NsamplesxNgridpoints.

###########
# Input required from continuousGaussianPlotting.R:
# D : sample points positions
# Resolution : the resolution of the spatial grid (high or low)
# sampleIDs2 : the list of samples to used (in this order)
# SD : the minimum sigma2 to use
###########

DIR="~/Documents/ClareDPhil/DPhil/Spain/densityMaps";#setwd(DIR)

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

# sample points and gridpoints in units of scaleFactor (10,000 metres)
X <- D[sampleIDs2,"Xmuni.ave.grand.precise"]
Y <- D[sampleIDs2,"Ymuni.ave.grand.precise"]

grand.points <- SpatialPoints(cbind(X,Y))
proj4string(grand.points) <- CRS(proj4string(espMapTotal)) # <-- coordinates of indiivuals come from original coordinates in espMapTotal
sp.grand.points <- spTransform(grand.points,CRS(proj4string(espMapTotalProj)))  # in grid coordinates

SP_grdPoints = SP_grd; gridded(SP_grdPoints) = FALSE  # convert grid to points (centres of grids)

scaleFactor=10000
samplePoints = sp.grand.points@coords/scaleFactor
#gridPoints = sp_grd2Points@coords/scaleFactor
gridPoints = SP_grdPoints@coords/scaleFactor


# Pick bandwidth such that sum of all sample points (weighted by a gaussian without variance-scaling) is some fixed value.
costFunction <- function(x,gridPoint=gridPoints[1,],sp=samplePoints,k=1){
  # sum over gaussian density at grid-point, with sigma2=x
  #sum(dmvnorm(sp,mean=gridPoint,sigma = x*diag(2)))
  #abs(sum(dmvnorm(sp,mean=gridPoint,sigma = x*diag(2))) - k)
  #sum(dmvnorm(sp,mean=gridPoint,sigma = x*diag(2))) - k
  #abs(sum(dmvnorm(sp,mean=gridPoint,sigma = 10*x*diag(2))) - k)  # factor of 10 helps with precision problems as sum of densities gets very small
  #abs(sum(dmvnorm(sp,mean=gridPoint,sigma = 100*x*diag(2))) - k)
  abs(sum( exp( -(1/(2*x))*(colSums((t(sp) - gridPoint)^2 )) ) ) - k)
}

# 61, 251,669,1026,16520,6213 ==> difficult grid position
# Test an example.

gp = cbind(gridPoints,startPoints)[6213,]
opt = optim(gp[3],costFunction,gridPoint=gp[1:2],sp=samplePoints,k=1)
X=seq(0,20,by=0.01) # values of sigma
plot(X,sapply(X,costFunction,gridPoint=gp[1:2],sp=samplePoints,k=1))
abline(v=opt$par)
abline(h=opt$value)

# Optimise for all gridpoints!
# pick a sensible start point based on an easy sum
dens = apply(gridPoints,1,function(gp){
  costFunction(SD,gridPoint=gp[1:2],sp=samplePoints,k=0)
})
  
# Pick k such that the highest density region has optimal bandwidth parameter SD
maxPoint = which(dens==max(dens))
gp=cbind(gridPoints,startPoints)[maxPoint,]
k = costFunction(SD,gridPoint=gp[1:2],sp=samplePoints,k=0) # k such
#h=sapply(,function(k) optim(gp[3],costFunction,gridPoint=gp[1:2],sp=samplePoints,k=k)$par)
#plot(0:50,h)
print(k)

# Find optimal bandwidth.
h = apply(cbind(gridPoints,startPoints),1,function(gp){
  print(date())
  opt= optim(gp[3],costFunction,gridPoint=gp[1:2],sp=samplePoints,k=k)
  return(c(opt$value,opt$par,opt$convergence))
})


# get list of gridpoints that are not in the sea
Notinthesea = which(!is.na(over(SP_grd,espMapTotalSansCanProj)[,1]))

newSDs = h[2,]
hist(newSDs[Notinthesea],breaks=100)
plot(startPoints,h[2,])

# Plot how the bandwidth changes!
sdPlot = plotDenistyContinuous(newSDs,colourSet=mapsColourSet,title="",sp_grd=SP_grd,topAlpha=0.5,fixedLims=c(0,max(newSDs[Notinthesea])),minDensity=0,
  transparencyVector=transparencyVector,notinthesea = Notinthesea,mainMap="SpainMap_under2Proj")

sdPlot[[2]] + spplot(sp.grand.points,col.regions="white",cex=0.4)


# Evaluate Gaussian density for each samplepoint/gridpoint pair with new bandwdith.
d = sapply(1:nrow(gridPoints),function(i){
    dmvnorm(samplePoints,mean=gridPoints[i,],sigma = (newSDs[i])*diag(2))
  })

d2 = sapply(1:nrow(gridPoints),function(i){
    t=dmvnorm(samplePoints,mean=gridPoints[i,],sigma = (newSDs[i])*diag(2))
    t*(2*pi*(newSDs[i]))
  })


##### Save the output!!
# h = variable bandwidth
# d = weights based on bandwidth for each grid point and sample point pair
# k = minimum bandwidth used.

save(d,h,k,SD,sampleIDs2,SP_grd,file=paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-weights-variableBWexpBalloon-',SD,'-',resolution,'Res.Rdata'))
# For weights without basque1
#save(d,h,k,SD,sampleIDs2,SP_grd,file=paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-weights2-variableBWexpBalloon-',SD,'-',resolution,'Res.Rdata'))



############### 
# check some points  

myPoints=c(6000:6213,13500:13600)
plot(1:length(newSDs[Notinthesea]),newSDs[Notinthesea],ylim=c(0,50))
plot(myPoints,newSDs[myPoints])

plot(myPoints,-log(colSums(dens)[myPoints]))
plot(myPoints,startPoints[myPoints])
plot(startPoints[myPoints],newSDs[myPoints])
abline(0,1)


for(g in myPoints[1:5]){
  gp =gridPoints[g,]
  opt = optim(10,costFunction,gridPoint=gp[1:2],sp=samplePoints,k=1)
  #opt = optim(1,costFunction,gridPoint=gp[1:2],sp=samplePoints,k=1,lower=0,method = "L-BFGS-B")
  X=seq(0,20,by=0.1) # values of sigma
  plot(X,sapply(X,costFunction,gridPoint=gp[1:2],sp=samplePoints,k=1))
  print(opt$par)
  abline(v=opt$par)
  abline(h=opt$value)
}

#myGridPoint = newSDs
myGridPoint = startPoints
myGridPoint[-myPoints] = 0
sdPlot = plotDenistyContinuous(myGridPoint,colourSet=mapsColourSet,title="",sp_grd=sp_grd2,topAlpha=0.5,fixedLims=NULL,minDensity=0,
  transparencyVector=transparencyVector,notinthesea = Notinthesea,mainMap="SpainMap_under2Proj")
  sdPlot[[2]] + spplot(sp.grand.points,col.regions="pink",cex=0.4) 
  

