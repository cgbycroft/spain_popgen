################# Plot mixture model - version 2 #####################

###### Global Vars ##########

setwd("~/Documents/ClareDPhil/DPhil/Spain/mixtureModel")
source("../Other/ClaresFunctions.R")
source("../mixtureModel/nnls.R")
colourSet <- c("white","red","blue","black")
some.colorsEnd2<-colorRampPalette(colourSet,interpolate="linear")(100) # as above, but with a dark grey final for capped values
pngHeight = 10000
G <- UnderMap2Proj + update(commapborders2Proj,col.regions=hsv(0,0,0.8),alpha.regions=0.5) + provmapborders2Proj + update(UnderMap2Proj,col.regions="transparent")
imageHeight=1500

# load tree structures for main spain analysis
load(paste0("../chromopainter/SPAIN.A2.ProcessedResults_v1-Spain_v1_v7a.new.RData"),verbose=TRUE)

###### Versions ##########

version_string = "Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9" # Spain splits from context fs results
version_string = "Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9-excHighBasFrance" # <==== This is what we used for the paper... Spain splits from context fs results
version_string = "Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9-excHighBasFrance-means" # 
version_string = "Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12" # <==== with Basques as donors. 
version_string = "Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine" # <==== with Basques as donors. For the paper, this was used for the GT analysis + mixture model!
version_string = "LOO-SP_PR_HM_NA_v10_DONORpopn.v12-excHighBasFrance" # <==== Donor groups only, including basques. made using mixture_model_donors.R. Used for Paper!!
###################

load(paste("../mixtureModel/mixtureModelOutput.",version_string,".Rdata",sep=""),verbose=T)

directory <- paste("/plots.",version_string,"_DIR",sep="")
system(paste('mkdir context_plots',directory,sep=""))

###### Update colours ######
recieverLegend2[c("Portugal","SpainPopn_100"),"colour"] = "darkblue"
recieverLegend2[c("Portugal","SpainPopn_100"),"shape"] = 22
recieverLegend2[c("Portugal","SpainPopn_100"),"Factor"] = "SpainPopn_100"
#recieverLegend2[c("France","donorPopn_17"),"colour"] = "yellow"
#recieverLegend2[c("NorthMorocco","donorPopn_92"),"colour"] = "blue"
#recieverLegend2[c("WesternSahara","donorPopn_90"),"colour"] = "red"
#recieverLegend2[c("donorPopn_80","Sub-saharan.Nigeria_80"),"colour"] = "orange"
recieverLegend2[c("donorPopn_318"),"colour"] = "cyan"
recieverLegend2[c("donorPopn_318"),"shape"] = 24
recieverLegend2[c("donorPopn_318"),"Factor"] = "donorPopn_318"
recieverLegend2[c("donorPopn_27"),"shape"] = 24
recieverLegend2[c("donorPopn_6"),"shape"] = 21
recieverLegend2[c("SpainPopn_100"),"shape"] = recieverLegend2[c("donorPopn_6"),"shape"] 
recieverLegend2[c("SpainPopn_100"),"colour"] = recieverLegend2[c("donorPopn_6"),"colour"] 

# Save this recieverLegend2
#save(recieverLegend2,file="../mixtureModel/myReceiverLegend2.RData")
load(file="../mixtureModel/myReceiverLegend2.RData",verbose=TRUE)

if(version_string=="Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine")  {
  recieverLegend2[paste0("SpainPopn_",colnames(numToShape)),"shape"] = numToShape[2,]
  recieverLegend2[paste0("SpainPopn_",colnames(numToShape)),"colour"] = numToColour[2,]
  numToColour = cbind(numToColour,c("100",recieverLegend2["SpainPopn_100","colour"])) # portuguese colours!
  colnames(numToColour)[ncol(numToColour)]="100"
}

###### Fix Netherlands typo ########## <== now fixed in getPopnLabels function.
#names(popnLabels)[names(popnLabels)=="Neatherlands-Germany"] = "Netherlands-Germany"

###### Other objects ##########

donorPopns <- colnames(mixmat)
SpainPopns <- grep("SpainPopn",rownames(mixmat),value=T)
ALLPopns <- rownames(SUMMARYmat)

  
##### Get data
data <- mixmat
dataS <- data[grepl("Spain",rownames(mixmat)),]
boot95=percentile95bootstrap
boot05=percentile05bootstrap

##### non-zero donors
nonZero1 = apply(boot05[SpainPopns,],2,function(x) min(x)>0) # greater than 0 in all SPain
nonZeroDonors1 = names(nonZero1)[nonZero1]
nonZero2 = apply(boot05[SpainPopns,],2,function(x) max(x)>0) # greater than 0 in at least one
nonZeroDonors2 = names(nonZero2)[nonZero2]
nonZero3 = apply(mixmat[SpainPopns,nonZero2],2,function(x) max(x)>0.001 ) # max estimate greater than 0.001 in at least one
nonZeroDonors3 = names(nonZero3)[nonZero3]
#nonZeroDonors4 = names(nonZero3&nonZero2)[nonZero3&nonZero2] # greater than zero in at least one, and no-zero estimate > 0.001

# total contribution
  
##### Set order of recipient populations
popOrder = names(sort(mixmat[grepl("Spain",rownames(mixmat)),"donorPopn_92"])) # sort by amount of nmoroccan
if( "donorPopn_318" %in% colnames(mixmat)){
  
  popOrder = names(sort(mixmat[grepl("Spain",rownames(mixmat)),"donorPopn_318"])) # sort by amount of basque  
  donorOrder = colnames(dataS)[order(apply(dataS[popOrder,],2,max))]
  # OR order by fineSTRUCTURE tree
  treeOrder = paste0("SPAIN_A2_",labels(tdendReorder))  
  #subTdendReorder = drop.leaves(tdendReorder,labels(tdendReorder)[!treeOrder%in%rownames(SPAIN.split.matrix)])
  clusterOrder = unique(SPAIN.split.matrix[treeOrder[treeOrder%in%rownames(SPAIN.split.matrix)],2]  )
  popOrder = paste0("SpainPopn_",c(100,clusterOrder))
  popOrder = popOrder[popOrder%in%rownames(mixmat)]
}
popOrder = rev(popOrder)

##### Set order of donor populations (based on donor trees...)
theseNA <- unlist(popnLabels[1:6])
theseEUR <- unlist(popnLabels[11:22])
theseSS <- unlist(popnLabels[7:10])

donorOrder=c(rev(rownames(mixmat)[rownames(mixmat)%in%paste0("donorPopn_",theseSS)]),
             rev(rownames(mixmat)[rownames(mixmat)%in%paste0("donorPopn_",theseNA)]),
          rev( rownames(mixmat)[rownames(mixmat)%in%paste0("donorPopn_",theseEUR)]))
donorOrder = c(donorOrder,donorPopns[!donorPopns%in%donorOrder])



###########################
##### basic bar plots
###########################

barPlotMatrix(data[c(popOrder,donorOrder),donorOrder],subset="post-mix-All",bootstraps=T,limit=0.005,dir=directory)
barPlotMatrix(data[donorOrder,donorOrder],subset="post-mix-Donors",bootstraps=T,limit=0.001,PDF=TRUE,
              dir=directory,width=1500,height=1500,colourCols=TRUE,oma=c(11,11,1,0),cex=0.4,yTextMargin = 1.25,labelStyle = 1)
barPlotMatrix(data[donorOrder,donorOrder],subset="post-mix-Donors",bootstraps=T,limit=0.001,PDF=FALSE,
              dir=directory,width=1500,height=1500,colourCols=TRUE,oma=c(11,11,1,0),cex=0.4,yTextMargin = 1.25,labelStyle = 1)

if("Portugal_6"%in%popOrder) barPlotMatrix(dataS[c("SpainPopn_49",popOrder[-1]),donorOrder],subset="post-mix-Spain",bootstraps=T,limit=0.005,dir=directory) else +
barPlotMatrix(dataS[popOrder,donorOrder],subset="post-mix-Spain",bootstraps=T,limit=0.001,dir=directory,width=1200,height=1200,yTextMargin=0.5,printLabs=F)
barPlotMatrix(dataS[popOrder,donorOrder[donorOrder!="donorPopn_17"]],subset="post-mix-Spain-exclFrance",bootstraps=T,limit=0.001,dir=directory,width=1200,height=1200,yTextMargin=0.5,xTextMargin=0.1,printLabs=F)
barPlotMatrix(dataS[popOrder,donorOrder[donorOrder%in%paste('donorPopn_',theseNA,sep="")]],subset="post-mix-Spain-NA",bootstraps=T,limit=0.001,dir=directory,width=1200,height=1200,yTextMargin=0.5,xTextMargin=0.1,printLabs=F)
barPlotMatrix(dataS[popOrder,donorOrder[donorOrder%in%nonZeroDonors1]],subset="post-mix-Spain-nonZeroAll",bootstraps=T,limit=0.001,dir=directory,width=1200,height=1200,yTextMargin=0.5,xTextMargin=0.1,printLabs=F)
barPlotMatrix(dataS[popOrder,donorOrder[donorOrder%in%nonZeroDonors2]],subset="post-mix-Spain-nonZeroAtLeastOne",bootstraps=T,limit=0.001,dir=directory,width=1200,height=1200,yTextMargin=0.5,xTextMargin=0.1,printLabs=F)


###########################
#### Gaussian maps!
###########################

stuff = getFilteredPointsSpain(filePrefix=version_string,
                               excludePont=T,split.matrix=SPAIN.split.matrix,Portugal=T,addSPAINA2=T)
D = stuff$D
D = D[rownames(D)%in%rownames(SPAIN.split.matrix),]

donors = unique(DONOR.split.matrix[,2])
summaryData = mixmat_inds[rownames(SPAIN.split.matrix),] # mixture model on individual level
colorSet = c("yellow","red","blue","black")
pointsPositions = D
split.matrix=SPAIN.split.matrix
symbols = recieverLegend2[paste0("SpainPopn_",SPAIN.split.matrix[,2]),"shape"]; names(symbols)=rownames(SPAIN.split.matrix)
addTransparency=FALSE
scaleLab="Fraction contribution"
#setSDbyData=FALSE
setSDbyData=TRUE
#model="t-dist"
SD=3.5
model="Gaussian"
resolution="high"

for(donor in donors){
  donorNames = getPopnLabel(paste0("donorPopn_",donor),popnLabels)
  print(donorNames[1])
  # input variables for gaussian plotting:
  inputVector = summaryData[,paste0("donorPopn_",donor)]
  names(inputVector) = rownames(summaryData)
  filename = paste0("context_plots",directory,"/contributionFromTEST-",donorNames[2])
  fixedLims = c(0,max(inputVector[sampleIDs[sampleIDs%in%names(inputVector)]]))
  title=donorNames[1]
  PlotPortugalToo=TRUE
  # actual plotting
  source('../densityMaps/continuousGaussianPlotting.R')
}

# plot raw copying means!
summaryData = SUMMARYmat_1mean[rownames(SPAIN.split.matrix),] # raw mean copying values for each individual
#setSDbyData=FALSE

for(donor in donors){
  donorNames = getPopnLabel(paste0("donorPopn_",donor),popnLabels)
  print(donorNames[1])
  # input variables for gaussian plotting:
  inputVector = summaryData[,paste0("donorPopn_",donor)]
  names(inputVector) = rownames(summaryData)
  filename = paste0("context_plots",directory,"/rawMeanCoancestryContributionFrom-",donorNames[2])
  fixedLims = c(min(inputVector),max(inputVector))
  title=donorNames[1]
  scaleLab="Mean coancestry (cM)"
  PlotPortugalToo=TRUE
  # actual plotting
  source('../densityMaps/continuousGaussianPlotting.R')
}

#############################################
# Spatially smoothed ancestry profiles.
#############################################
#filter = 80

# Load weights matrix d from getVariableBandwidth.R. NOTE: weights depend on the set of sampled points used. i.e sampleIDs2.
if(version_string=="Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine") {
   # weights based on 676 spanish individuals without Basque1
  load(paste0('../densityMaps/Spain.A2-SpatialDistances-',filter,'Km-weights2-variableBWexpBalloon-',SD,'-',resolution,'Res.Rdata'),verbose=TRUE)
  
  } else {
    # based on 726 spanish individuals
    load(paste0('../densityMaps/Spain.A2-SpatialDistances-',filter,'Km-weights-variableBWexpBalloon-',SD,'-',resolution,'Res.Rdata'),verbose=TRUE)
  }

# SUMMARYmat is group-level matrix to be populated
# SUMMARYmat_1 is a matrix of coancestry-vectors (chunklengths summed within each donor group). size nxK where n is number of individual samples, and K is number of donor groups

sumd <- colSums(d)
newmat <- SUMMARYmat[donorPopns,]
tmpmat <- newmat
#SUMMARYmat_inds <- rbind(SUMMARYmat_1[grep("SPAIN_A2",rownames(SUMMARYmat_1)),],newmat)
SUMMARYmat_locs <- t(t(SUMMARYmat_1[sampleIDs2,])%*%d)/sumd  # IMPORTANT: The distances matrix d is calculated as in the getVariableBandwidth.R script above and samples have order sampleIDs2
toFit = SUMMARYmat_locs/rowSums(SUMMARYmat_locs)  # Normalise coancestry vectors to sum to one.
                                
mixmat_locs <- matrix(0,nrow=nrow(SUMMARYmat_locs),ncol=nrow(newmat))
  colnames(mixmat_locs) <- donorPopns
  rownames(mixmat_locs) <- rownames(SUMMARYmat_locs)
mixmat_output_locs <- list()

# Each 'pop' here is a grid-point in Spain!
for(pop in 1:nrow(mixmat_locs)){
  if( pop%%1000==0 ) print(pop)
  ourmixALL <- getoverallfit(tmpmat/rowSums(tmpmat),
                        toFit[pop,])
  #mixmat_output_locs[[pop]] <- ourmixALL
  ourmix <- ourmixALL$x
  ourmix <- ourmix[ourmix>0]
 # ourmix <- ourmix/sum(ourmix)
  mixmat_locs[pop,colnames(mixmat_locs)%in%names(ourmix)] <- ourmix
}

# portugal - just one value (the mean)!
spainPoints = pointsPositions[rownames(pointsPositions)%in%rownames(SPAIN.split.matrix),c("Xmuni.ave.grand.precise","Ymuni.ave.grand.precise")]
    portPoints = pointsPositions[rownames(pointsPositions)%in%rownames(SPAIN.split.matrix),c("Xmuni.mixture","Ymuni.mixture")]
    colnames(portPoints) = colnames(spainPoints)
    myPoints = rbind(spainPoints[!is.na(spainPoints[,1]),],portPoints[grepl("POPRES",rownames(portPoints)),])
    inPort = over(makeSpPoints(myPoints),PortugalMap);
    inPort = rownames(inPort)[!is.na(inPort[,1])]
    inPort = inPort[inPort%in%rownames(SPAIN.split.matrix)]
    # treat each point as having equal weight ==> just the mean value!
    SUMMARYmat_locsPort = t(t(SUMMARYmat_1[inPort,])%*%cbind(rep(1,length(inPort)),rep(1,length(inPort))))/length(inPort)
toFitPort = SUMMARYmat_locsPort[1,]/rowSums(SUMMARYmat_locsPort)[1]
mixmat_locsPort = mixmat_locs[1,]
ourmixALL <- getoverallfit(tmpmat/rowSums(tmpmat),
                        toFitPort)
  #mixmat_output_locs[[pop]] <- ourmixALL
  ourmix <- ourmixALL$x
  ourmix <- ourmix[ourmix>0]
 # ourmix <- ourmix/sum(ourmix)
  mixmat_locsPort[names(mixmat_locsPort)%in%names(ourmix)] <- ourmix
  
#save(mixmat_locs,mixmat_locsPort,model,SD,file=paste0(version_string,"-",model,"-",SD,"-geolocalMixtureModel-",resolution,"Res.RData"))
#save(mixmat_locs,mixmat_locsPort,model,SD,file=paste0(version_string,"-",model,"-",SD,"-geolocalMixtureModel_variableBW-",resolution,"Res.RData"))
#save(mixmat_locs,mixmat_locsPort,model,SD,file=paste0(version_string,"-",model,"-",SD,"-geolocalMixtureModel_variableBWsamplebased-",resolution,"Res.RData"))
save(mixmat_locs,mixmat_locsPort,model,SD,file=paste0(version_string,"-",model,"-",SD,"-geolocalMixtureModel_variableBWexpBalloon-",resolution,"Res.RData"))  # <===== USED THIS ONE IN THE END!

###################################################### 
# plot this!

# some objects that would normally be created if had run continuousGaussianPlotting.R
load(paste0('~/Documents/ClareDPhil/DPhil/Spain/densityMaps/Spain.A2-SpatialDistances-',filter,'Km-highRes.Rdata'),verbose=T)  # this is the new projected version
Notinthesea = which(!is.na(over(SP_grd,espMapTotalSansCanProj)[,1]))
transparencyVector=NULL

#####
#bwmethod="variableBW"  
bwmethod="variableBWexpBalloon"
plainPoints=TRUE
#plainPoints=FALSE
#####

load(file=paste0(version_string,"-",model,"-",SD,"-geolocalMixtureModel_",bwmethod,"-",resolution,"Res.RData"),verbose=TRUE)

#mapsColourSet=viridis
mapsColourSet=magma

if(resolution=="high") myGrid = sp_grdZoom2 else  myGrid = sp_grd2

symbols = recieverLegend2[paste0("SpainPopn_",SPAIN.split.matrix[,2]),"shape"]; names(symbols)=rownames(SPAIN.split.matrix)
ptsize=0.8
pointsBorderCols = rep("black",length(symbols))
names(pointsBorderCols) = names(symbols)

if(version_string=="Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine"){
  # put grey borders on dark colours!
  bgCols = recieverLegend2[paste0("SpainPopn_",SPAIN.split.matrix[,2]),"colour"]
  dark = rgb2hsv(col2rgb(bgCols))["v",]<0.1
  pointsBorderCols[dark] = "gray"
  names(pointsBorderCols) = names(symbols)
}

if(plainPoints) {
  symbols[]=4
  #ptsize=0.5
  ptsize=0.7
  pointsBorderCols=NULL
}


for(donor in donors){  
  donorNames = getPopnLabel(paste0("donorPopn_",donor),popnLabels)
  print(donorNames[1])
  # input variables for gaussian plotting:
  inputVector = summaryData[,paste0("donorPopn_",donor)]
  names(inputVector) = rownames(summaryData)
  filename = paste0("context_plots",directory,"/GeoLocalCoancestryContributionFrom-",donorNames[2])
  title=donorNames[1]
  scaleLab="Fraction contribution"
  toPlot = mixmat_locs[,paste0("donorPopn_",donor)]
  toPlotPort = mixmat_locsPort[paste0("donorPopn_",donor)]
  
  if(sum(c(toPlot[Notinthesea],toPlotPort)>0)<2) next
  
  #fixedLims = c(min(c(toPlot,toPlotPort)),max(c(toPlot,toPlotPort))) # That includes gridpoints in the sea...
  myValues = c(toPlot[Notinthesea],toPlotPort)
  fixedLims = c(min(myValues[myValues!=0]),max(c(toPlot[Notinthesea],toPlotPort))) # Don't include zero in the fixed limits.
  fixedLims = c(0,max(c(toPlot[Notinthesea],toPlotPort))) # include zero (v3)
  
  Plots <- plotDenistyContinuous(toPlot,colourSet=mapsColourSet,title="",sp_grd=myGrid,topAlpha=0.5,fixedLims=fixedLims,transparencyVector=transparencyVector,notinthesea = Notinthesea,mainMap="SpainMap_under2Proj")
  PlotsPort <- plotDenistyContinuous(toPlotPort,colourSet=mapsColourSet,title="",sp_grd=myGrid,topAlpha=0.5,fixedLims=fixedLims,transparencyVector=transparencyVector,notinthesea = InPortOnly,mainMap="SpainMap_under2Proj")
  
  #sdPlot = plotDenistyContinuous(newSDs,colourSet=mapsColourSet,title="",sp_grd=myGrid,topAlpha=0.5,fixedLims=range(newSDs),transparencyVector=transparencyVector,notinthesea = Notinthesea,mainMap="SpainMap_under2Proj")

  # withPortugal
  G = update(SpainMap_under2Proj,col.regions="white") + Plots[[2]] + PlotsPort[[2]]  + update(commapborders2Proj,lwd=0.5,col=hsv(0,0,0.6))
  File = paste0(filename,"-",model,"-",SD,"-",resolution,"res-withPort.png")
  if(plainPoints) File = paste0(filename,"-",model,"-",SD,"-",resolution,"res-withPort-plainPoints3.png") # <= v3 means bigger crosses

  vectorForHist = as.vector(toPlot[Notinthesea])  # Histogram shows all values of the notinthesee grid-points (excluding portugal).
  vectorForHist = vectorForHist[vectorForHist>0] # only show histogram for non-zero values.
  png(File,width=imageHeight*spainAspect*2,height=imageHeight*2,res=2*150)
  par(cex=1.5,cex.lab=0.8)  # These affect the histogram
  plotPointsOnMaps2(D=pointsPositions[rownames(pointsPositions)%in%names(inputVector),],dend=NULL,zoom=F,split.matrix=SPAIN.split.matrix,m=2,
                              showLines=F,plotDendro=F,dendroPlot=NULL,baseMap=G,
                              levelsVector=NULL,shapesVector=symbols,levelSplits=10,
                              extraTitle="",pointSize=ptsize,LWD=2,minSplitSize=0,
                              title=title,vectorForHist=vectorForHist,
                              trans=1,darken=1,showScale=T,colourSet=mapsColourSet,
                              fixedLims=fixedLims,histXlims=fixedLims,extraPoints=NULL,
                              mapLabels=T,newPage=T,numToColours=numToColour,pointLWD=0.7,
                              histBgCol="transparent",jitter="max",morePoints=NULL,scaleLab=scaleLab,
                              histBinWidth=NULL,bbox=rbind(c(-9.5,-5),c(40,44)),pointsPCH=21,
                              pointsBorderCol=pointsBorderCols,defaultPointCol = "gray")
  dev.off()
  
  # without Portugal
  #fixedLims = c(min(c(toPlot[Notinthesea])),max(c(toPlot[Notinthesea])))
  G = update(SpainMap_under2Proj,col.regions="white") + Plots[[2]] + update(commapborders2Proj,lwd=0.5,col=hsv(0,0,0.6))
  File = paste0(filename,"-",model,"-",SD,"-",resolution,"res.png")
  if(plainPoints) File = paste0(filename,"-",model,"-",SD,"-",resolution,"res-plainPoints.png")
  
  vectorForHist = as.vector(toPlot[Notinthesea])
  vectorForHist = vectorForHist[vectorForHist>0] # only show histogram for non-zero values.
  png(File,width=imageHeight*spainAspect*2,height=imageHeight*2,res=150*2)
  par(cex=1.5,cex.lab=0.8)  # These affect the histogram
  plotPointsOnMaps2(D=pointsPositions[rownames(pointsPositions)%in%names(inputVector),],dend=NULL,zoom=F,split.matrix=SPAIN.split.matrix,m=2,
                              showLines=F,plotDendro=F,dendroPlot=NULL,baseMap=G,
                              levelsVector=NULL,shapesVector=symbols,levelSplits=10,
                              extraTitle="",pointSize=ptsize,LWD=2,minSplitSize=0,
                              title=title,vectorForHist=vectorForHist,
                              trans=1,darken=1,showScale=T,colourSet=mapsColourSet,
                              fixedLims=fixedLims,histXlims=fixedLims,extraPoints=NULL,
                              mapLabels=T,newPage=T,numToColours=numToColour,pointLWD=0.6,
                              histBgCol="transparent",jitter="max",morePoints=NULL,scaleLab=scaleLab,
                              histBinWidth=NULL,bbox=rbind(c(-9.5,-5),c(40,44)),pointsPCH=21,
                              pointsBorderCol=pointsBorderCols,defaultPointCol = "gray")
  dev.off()
  
}

##################### 
# plot main 4 components for each cluster, with variable heights
#############################################

    nTop = 10
    #s = colSums(dataS[,colnames(dataS)%in%c(nonZeroDonors1,nonZeroDonors2)])
    s = colSums(dataS[,colnames(dataS)%in%c(nonZeroDonors3)])
    ord = order(s,decreasing=T)
    i= ord
    #i = ord[1:nTop]
    #i=9
    #cols = c("#1E1175","#3D4ECC","#EE1505","#33663E","#EE1505","#EE1505")[1:4]
    cols = recieverLegend2[names(s)[i],"colour"]
    names(cols) = names(s)[i]
    
    #png(paste0('context_plots',directory,'/',version_string,"-barplot-top",nTop,"Donors.png"),
    #png(paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonorsTEST.png"),
    
    #png(paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonors.png"), 
    #bg="transparent",width=950,height=1300,res=150)
    pdf(paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonors.pdf"), 
    bg="transparent",width=950/150,height=1300/150)
    
    par(mfrow=c(length(i)+1,1),mar=c(1,1,1,1),mgp=c(1,1,0),oma=c(7,14,1,1))
  
    #png(paste0('context_plots',directory,'/',version_string,"-barplot-Basque.png"),     
    #bg="transparent",width=1600,height=1000,res=150)
    #par(mfrow=c(length(i)+1,1),mar=c(1,5,1,1),mgp=c(1,1,0))
    

    nam = sapply(names(cols),getPopnLabel,popnLabels)[1,]
  
    for(donor in names(s)[i]){
          a = mixmat[popOrder,donor]
          b = barplot(a,names.arg=NA,col=cols[donor],ylim=c(0,max(boot95[popOrder,donor])),
                      cex.lab=1,las=1,yaxt="n")
          text(-0.7,max(boot95[popOrder,donor])/2,nam[donor],xpd=NA,pos=2,cex=2)
          axis(2,at = c(0,max(boot95[popOrder,donor])),cex=2,labels=c("",round(max(boot95[popOrder,donor]),3)),las=2)
          #axis(2,at=max(boot95[popOrder,donor])/2,labels=round(max(boot95[popOrder,donor]),5),las=2,tick=FALSE)
          segments(b,boot05[popOrder,donor],b,boot95[popOrder,donor],xpd=NA)
          segments(b-0.1,boot05[popOrder,donor],b+0.1,boot05[popOrder,donor],xpd=NA)
          segments(b-0.1,boot95[popOrder,donor],b+0.1,boot95[popOrder,donor],xpd=NA)            
          
         points(x=-0.5,y=par()$usr[4]/2,xpd=NA,
         bg=recieverLegend2[donor,"colour"],
         col="black",
         pch=recieverLegend2[donor,"shape"],
         xpd=NA,cex=2)

    }
  xline = strheight("M", units = "user")
  #print(xline)
  #abline(h = xline)
  mtext("Iberian clusters",xpd=NA,side=1,line=4.5,cex=1.2,col=add.alpha("black",0.5),outer = TRUE)
  mtext("Donor groups",xpd=NA,side=2,line=8,outer = TRUE,cex=1.2,col=add.alpha("black",0.5))
  
  points(b,y=rep(par()$usr[3] - 3*xline,length(b)),
         bg=recieverLegend2[popOrder,"colour"],
         col="black",
         pch=recieverLegend2[popOrder,"shape"],
         xpd=NA,cex=4)

  # Add labels
  nGroups=length(popOrder)
  #mtext("Iberian clusters",side=1,line=3.5,cex=1.2,col=add.alpha("black",0.5))
  xline = strheight("M", units = "user")
  par(srt=45)
  text(x=b,y=rep(par()$usr[3] - 5*xline,nGroups),cex=2,pos=2,xpd=NA,labels=sapply(popOrder,getPopnLabel,popnLabels)[1,])
 
  if("SpainPopn_100"%in%popOrder) {
      mtext("Portugal",side=1,line=3,cex=1,las=2,at = b[which(popOrder=="SpainPopn_100")])
  
  }
    
  dev.off()


  t=mixmat[popOrder,names(s)[ord]]
  nam = sapply(names(s)[ord],getPopnLabel,popnLabels)[1,]
  colnames(t) = nam
  rownames(t) = sapply(rownames(t),getPopnLabel,popnLabels)[1,]
    
  
###########################
#### Plot mixture model as pies.
###########################
    nTop = 10
    s = colSums(dataS[,colnames(dataS)%in%c(nonZeroDonors3)])
    ord = order(s,decreasing=T)
    i= ord
    nam = sapply(names(cols),getPopnLabel,popnLabels)[1,]
  
# plot locations as pies
# get centres for pies

centreOrder = sapply(popOrder,getPopnLabel,popnLabels)[1,]
pieCentres = as.data.frame(rbind(c(3967776,3452452),c(3911320,2963600),c(4243934,2863600),c(4311269,3212974),c(4773060,3351395),c(4411320,3484428)))
colnames(pieCentres) = c("X","Y")
rownames(pieCentres) = popOrder
pieCols = recieverLegend2[names(s)[i],"colour"]; names(pieCols)=names(s)[i]

dummyVector = rep(0,length(pieCols)); names(dummyVector) = names(s)[i]

pieProps = t(mixmat[popOrder,names(s)[i]])
imageHeight=10
baseMap = UnderMap2Proj #+ update(commapborders2Proj,col.regions=hsv(0,0,0.8),alpha.regions=0.5) + update(UnderMap2Proj,col.regions="transparent")

File = paste0('context_plots',directory,'/',version_string,'-mixmatNonZeroDonors-pies.png')
pdf(gsub(".png",".pdf",File),height=imageHeight*contextAspect,width=imageHeight,bg="transparent")
plot.new()
p = piesOnMaps2(props=pieProps,pieCentres=pieCentres,baseMap=baseMap,pieColours=pieCols,
                     radius=7e4,adj=c(1.3,1),dummyMap=UnderMap2Proj,mapLabels=NULL)
p
dev.off()





###########################
### Print a table for latex

    s = colSums(dataS[,colnames(dataS)%in%c(nonZeroDonors3)])
    ord = order(s,decreasing=T)
    
donorsHere = names(s)[ord]
pointEst  = mixmat[popOrder,donorsHere]
b05 = boot05[popOrder,donorsHere]
b95 = boot95[popOrder,donorsHere]

forPrint = round(pointEst,4)
for(i in donorsHere){
  forPrint[,i] = paste0(forPrint[,i],"  (",round(b05[,i],4)," - ",round(b95[,i],4),")")
}

# Make thumbnail images of points
pdf(paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonors-DONOR-THUMBS.pdf"),height=1,width=1,bg="transparent")
par(mar=c(0,0,0,0),mgp=c(0,0,0))
for(i in 1:length(donorsHere)){
  print(i)
  pop = donorsHere[i] 
  stuff = recieverLegend3[gsub(".null","",pop),]
  plot.new()
  points(0.5,0.5,col="black",lwd=1.5,bg=stuff[["colour"]],pch=stuff[["shape"]],cex=7)
}
dev.off()

pdf(paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonors-RECIPIENT-THUMBS.pdf"),height=1,width=1,bg="transparent")
par(mar=c(0,0,0,0),mgp=c(0,0,0))
for(i in 1:length(popOrder)){
  print(i)
  pop = popOrder[i] 
  stuff = recieverLegend3[gsub(".null","",pop),]
  plot.new()
  points(0.5,0.5,col="black",lwd=1.5,bg=stuff[["colour"]],pch=stuff[["shape"]],cex=7)
}
dev.off()

# add in image files to latex table
for(i in 1:length(donorsHere)){
  colnames(forPrint)[i] = paste0("\\includegraphics[width=\\baselineskip,height=\\baselineskip,page=",i,",valign=t]{{{/Users/clare/Documents/ClareDPhil/DPhil/Spain/mixtureModel/context_plots",directory,"/",version_string,"-barplot-nonZeroDonors-DONOR-THUMBS}}} ",getPopnLabel(colnames(forPrint)[i],popnLabels)[1])
}
for(i in 1:length(popOrder)){
  rownames(forPrint)[i] =  paste0("\\includegraphics[width=\\baselineskip,height=\\baselineskip,page=",i,",valign=t]{{{/Users/clare/Documents/ClareDPhil/DPhil/Spain/mixtureModel/context_plots",directory,"/",version_string,"-barplot-nonZeroDonors-RECIPIENT-THUMBS}}} ",getPopnLabel(rownames(forPrint)[i],popnLabels)[1])

}
  
# Rename some columnts
#colnames(forPrint)[1:11]= c("Target population","Sample size","One date (generations)","One date (year)","Two dates – date 1 (generations)","Two dates – date 1 (year)","Two dates – date 2 (generations)","Two dates – date 2 (year)","fit.quality.1event","maxR2fit.1date",	"maxScore.2events")
write.table(forPrint,sep=",",file=paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonors.txt"),quote=FALSE,col.names=TRUE,row.names=TRUE)


S.sizes = table(paste0("SpainPopn_",SPAIN.split.matrix[,2]))
S.sizes = S.sizes[popOrder]

# For latex
xtab = xtable(t(forPrint))
# ===> Manually check widths in latex!
#align(xtab) = myAlign(forPrint[,1:11],widths=c(0.1,0.045,0.06,0.095,0.06,0.095,0.06,0.095,0.082,0.082,0.082),totalWidth = "\\linewidth")
digits(xtab) = xdigits(xtab)
# Make extra row for extra headings

sink(file = paste0('context_plots',directory,'/',version_string,"-barplot-nonZeroDonors.tex"))
print(xtab,size="footnotesize",floating=FALSE,
      include.rownames = TRUE,sanitize.colnames.function = bold,booktabs = FALSE)
sink()

  
################ Plot a simple dendrogram to link all the clusters
tdendReorderPhy = dend2phylo(tdendReorder)
te = extract.clade(tdendReorderPhy,node=1414)
trick = SPAIN.split.matrix;colnames(trick)[1] = "Split_49"
if(sum(grepl("SPAIN",rownames(trick)))>0 ) rownames(trick) = gsub("SPAIN_A2_","",rownames(trick))
hmm = getColouredDendro(tdendReorderPhy,split.matrix=trick,m=27,collapse=T,colLevel=49,returnDend=F)


##################### 
# plot residuals
##################### 

residuals <- vapply(mixmat_output[grep("Spain",names(mixmat_output))],FUN=function(x) residuals(x)[,1],FUN.VALUE=vector('double',length=length(grep("donor",rownames(mixmat)))))
hist(residuals,breaks=50)
# NOTE: residuals are observed - fitted

# compute Total Variation Distance between actual and modelled vectors

TVD=function(a,b){
  0.5*sum(abs(a-b))
}

TVDmix <- sapply(grep("Spain",names(mixmat_output),value=T),FUN=function(i) {
  y=mixmat_output[[i]]
  print(i) 
  coef=y$x
  dons = SUMMARYmat[names(y$x),names(y$x)]
  a = (t(dons)%*%coef)[,1]; a = a/sum(a)
  obs=SUMMARYmat[i,colnames(SUMMARYmat)!=i]
  b=obs/sum(obs)
  TVD(a,b)
})

# raw residuals
File=paste0('context_plots/',directory,'/',version_string,'-Residuals.png')
toPlot <- t(residuals[rev(1:nrow(residuals)),])
#colours <- c(colorScale(seq(min(toPlot),0,length.out=20),colourSet=c("yellow","red","white"),nBreaks=20),
#             colorScale(seq(0,max(toPlot),length.out=20),colourSet=c("white","blue","black"),nBreaks=20))
plotMixtureHeat(toPlot=toPlot[popOrder,donorOrder],filename=File,title='mixture model residuals: inferred - actual',cap=F,
                height=0.5*pngHeight*ncol(toPlot)/nrow(toPlot),width=pngHeight,res=150,Cex=3)

# Line plot <=== used in supp. info
File=paste0('context_plots/',directory,'/',version_string,'-Residuals2.png')
toPlot2 = toPlot[popOrder,donorOrder]
#png(File,height=1000,width=1200,res=150)
pdf(gsub(".png",".pdf",File),height=1000/150,width=1200/150,bg="transparent")
par(mar=c(5,5,2,11),cex.lab=1.2)
plot(NULL,ylim=c(-max(toPlot2),max(toPlot2)),xlim=c(1,length(popOrder)),xaxt="n",xlab="Iberian clusters",ylab="Linear mixture residual")
for( i  in donorPopns){
  lines(1:length(popOrder),toPlot2[,i],col=recieverLegend2[i,"colour"],lwd=0.5,lty=3)
  points(1:length(popOrder),toPlot2[,i],bg=recieverLegend2[i,"colour"],
         col="black",pch=recieverLegend2[i,"shape"],cex=1.3,lwd=0.5)
}
#axis(1,at=1:length(popOrder),labels=sapply(popOrder,getPopnLabel,popnLabels)[3,],las=2)
points(1:length(popOrder),y=rep(par()$usr[3]-diff(par()$usr)[3]/30,length(popOrder)),cex=3,
       bg=recieverLegend2[popOrder,"colour"],col="black",pch=recieverLegend2[popOrder,"shape"],xpd=NA)
abline(h=0,lty=3)
donorOrderLegend = colnames(toPlot2)[order(toPlot2[popOrder[1],],decreasing = TRUE)]
legend("topright",title="Donor groups",inset=c(-0.45,-0.05),xpd=NA,legend=sapply(donorOrderLegend,getPopnLabel,popnLabels)[1,],bty="n",col="black",pt.bg=recieverLegend2[donorOrderLegend,"colour"],pch=recieverLegend2[donorOrderLegend,"shape"])
#text("France")
dev.off()


##################### 
# plot locations of donor populations
##################### 

withPort=TRUE
withPort=FALSE
withSpainAndPort=FALSE
#withSpainAndPort=TRUE

DONOR.split.matrix2 = DONOR.split.matrix

# add in Portuguese group and/or Spain
if(withPort|withSpainAndPort){
  # add Portuguese to the DONOR.split.matrix
  port = SPAIN.split.matrix[grep("POPRES",rownames(SPAIN.split.matrix)),]
  port[,2]=6; colnames(port)=colnames(DONOR.split.matrix)
  DONOR.split.matrix2  = rbind(DONOR.split.matrix,port)
  donorPopns = c(donorPopns,"donorPopn_6")
  donorOrder = c(donorOrder[1:which(donorOrder == "donorPopn_3")],"donorPopn_6",donorOrder[(which(donorOrder == "donorPopn_3")+1):length(donorOrder)])
  popnLabels$Portugal = c(popnLabels$Portugal,6)

  if(withSpainAndPort){
  spain = SPAIN.split.matrix[grep("SPAIN",rownames(SPAIN.split.matrix)),]
  spain[,2]=1000; colnames(spain)=colnames(DONOR.split.matrix)
  DONOR.split.matrix2  = rbind(DONOR.split.matrix2,spain)
  #donorPopns = c(donorPopns,paste0("SpainPopn_",unique(spain[,2])))
  }
} else {
  donorPopns = donorPopns[donorPopns!="donorPopn_6"]
  donorOrder = donorOrder[donorOrder!="donorPopn_6"]
}

splitColours <- rbind(gsub("donorPopn_|SpainPopn_","",recieverLegend2[donorPopns,"Factor"]),recieverLegend2[donorPopns,"colour"])
    colnames(splitColours)=splitColours[1,]
splitShapes <- rbind(as.numeric(gsub("donorPopn_|SpainPopn_","",recieverLegend2[donorPopns,"Factor"])),recieverLegend2[donorPopns,"shape"])
    colnames(splitShapes)=as.character(splitShapes[1,])
if(withSpainAndPort){
  splitColours = cbind(splitColours,c(1000,"black"))
  splitColours[2,"6"] = "black" # Make portuguese black too
  splitShapes = cbind(splitShapes,c(1000,16))
  colnames(splitShapes)[ncol(splitShapes)] = "1000" 
  colnames(splitColours)[ncol(splitColours)] = "1000"
  donorOrder = c(donorOrder,"SpainPopn_1000")
}
  
  D <- point.concordance.PR_HM_NA[which(rownames(point.concordance.PR_HM_NA)%in%rownames(DONOR.split.matrix2)),]
  stuff = getFilteredPointsSpain(filePrefix=version_string,split.matrix = DONOR.split.matrix2,
                                excludePont=T,Portugal=TRUE,addSPAINA2=T)
  extras = stuff$D[(rownames(stuff$D)%in%c(rownames(DONOR.split.matrix2)))&(!rownames(stuff$D)%in%rownames(D)),colnames(stuff$D)%in%colnames(D)]
  # Exclude the portuguese because they're not randonly there
  if(withSpainAndPort) D = D[!rownames(D)%in%rownames(extras),]
  
  if(nrow(extras)>0) {
  #  rownames(extras) <- paste("SPAIN_A2_",rownames(extras),sep="")
  D2 <- rbind(D[,colnames(extras)],extras)
    
  } else {D2 = D}

  
# remove lone Portuguese individual!???? (nah, don't do that.)
#if(!withPort) D2 = D2[rownames(D2)!="POPRES_47614",]
imageHeight=10


# remove unknown location of Basque individuals
#noLocBasq = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==318][D2[rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==318],3]==0]
#D2 = D2[!rownames(D2)%in%noLocBasq,]
#D2 = D2[rownames(D2)!="SPAIN_A2_16176",] # remove spanish guy outside the main map

# Sensible plotting order!
OrderPlot = order.by.number.occurrences(DONOR.split.matrix2[rownames(D2),2])

File = paste0('context_plots',directory,'/',version_string,'-DonorLocations.png')
if(withPort) File = gsub('.png','-withPort.png',File)
if(withSpainAndPort) File = gsub('.png','-withPortandSpain.png',File)

#png(File,height=imageHeight*450*contextAspect,width=imageHeight*450,res=450)
pdf(gsub(".png",".pdf",File),height=imageHeight*contextAspect,width=imageHeight)
plotPointsOnMaps2(D2[OrderPlot,],split.matrix=DONOR.split.matrix2,m=2,baseMap=ALL_context_map,pointSize=1,LWD=1,
                              minSplitSize=15,
                              title=paste(length(donorPopns)," groups"),
                              mapLabels=F,newPage=T,
                              numToColours=splitColours,
                              numToShapes=splitShapes,pointLWD=0.1)
dev.off()

# All dots black
splitColours2 = splitColours; splitColours2[2,]=My.add.alpha("black",0.5)
splitShapes2 = splitShapes; splitShapes2[2,]=16

pdf(gsub(".png","BLACK.pdf",File),height=imageHeight*contextAspect,width=imageHeight)
plotPointsOnMaps2(D2[OrderPlot,],split.matrix=DONOR.split.matrix2,m=2,baseMap=ALL_context_map,pointSize=1,LWD=1,
                              minSplitSize=15,
                              title=paste(length(donorPopns)," groups"),
                              mapLabels=F,newPage=T,
                              numToColours=splitColours2,
                              numToShapes=splitShapes2,pointLWD=0.1)
dev.off()

#### Plot each group separately with others washed out.
File = gsub('.png','-groupsSeparate.png',File)
OrderPlot = order.by.number.occurrences(DONOR.split.matrix2[rownames(D2),2])
  
pdf(gsub(".png",".pdf",File),height=imageHeight*contextAspect,width=imageHeight)
for(donor in c(318,92)){  
  donorNames = getPopnLabel(paste0("donorPopn_",donor),popnLabels)
  print(donorNames[1])
  
  splitColours2 = splitColours; splitColours2[2,colnames(splitColours2)!=as.character(donor)]=My.add.alpha(splitColours2[2,colnames(splitColours2)!=as.character(donor)],0.1)
  splitShapes2 = splitShapes; 
  These = which(DONOR.split.matrix2[rownames(D2),2]==donor)
  OrderPlot2 = c(These,OrderPlot[!OrderPlot%in%These])
  
  plotPointsOnMaps2(D2[OrderPlot2,],split.matrix=DONOR.split.matrix2,m=2,baseMap=ALL_context_map,pointSize=1,LWD=1,
                                minSplitSize=15,
                                title=paste(length(donorPopns)," groups"),
                                mapLabels=F,newPage=T,
                                numToColours=splitColours2,
                                numToShapes=splitShapes2,pointLWD=0.1)
}
dev.off()

#### Plot Basque and North morocco 
File = gsub('.png','-BasqueNorthMorocco.png',File)

pdf(gsub(".png",".pdf",File),height=imageHeight*contextAspect,width=imageHeight)
  others = 
  splitColours2 = splitColours; splitColours2[2,!colnames(splitColours2)%in%c("318","92")]=My.add.alpha(splitColours2[2,!colnames(splitColours2)%in%c("318","92")],0.2)
  splitShapes2 = splitShapes; splitShapes2[2,!colnames(splitShapes2)%in%c("318","92")]=NA

  OrderPlot = order.by.number.occurrences(DONOR.split.matrix2[rownames(D2),2])

  plotPointsOnMaps2(D2[OrderPlot,],split.matrix=DONOR.split.matrix2,m=2,baseMap=ALL_context_map,pointSize=0.8,LWD=1,
                                minSplitSize=15,
                                title=paste(length(donorPopns)," groups"),
                                mapLabels=F,newPage=T,
                                numToColours=splitColours2,
                                numToShapes=splitShapes2,pointLWD=0.1)
dev.off()


# plot locations as pies
# pick groups
pops <- donorPopns

# get centres for pies
STUFF = list("D"=D2)
countries = unique(D2$geocode_country)
#pieCentres = getPopnCentres(donorPopns,DONOR.split.matrix,STUFF,columns=c("Xcountry","Ycountry"))
pieCentres = as.data.frame(t(sapply(countries,function(s) c(D2$Xcountry[D2$geocode_country==s][1],D2$Ycountry[D2$geocode_country==s][1]))),stringAsFactors=FALSE)
colnames(pieCentres) = c("X","Y")
pieCols = recieverLegend2[donorPopns,"colour"]; names(pieCols)=donorPopns
dummyVector = rep(0,length(donorPopns)); names(dummyVector) = donorPopns
pieProps = sapply(countries,function(s){
  t=table(DONOR.split.matrix[rownames(D2)[D2$geocode_country==s],2])
  p = t/sum(t)
  p2 = dummyVector; p2[paste0("donorPopn_",names(p))]=p;
  return(as.vector(p2))
})
rownames(pieProps) = names(dummyVector)

File = paste0('context_plots',directory,'/',version_string,'-DonorLocations-pies.png')
if(withPort) File = gsub('.png','-withPort.png',File)
#png(File,height=imageHeight*150*contextAspect,width=imageHeight*150,res=150)
pdf(gsub(".png",".pdf",File),height=imageHeight*contextAspect,width=imageHeight,bg="transparent")
plot.new()
p = piesOnMaps2(props=pieProps,pieCentres=pieCentres,baseMap=ALL_context_map,pieColours=pieCols,
                     radius=1,adj=c(1.3,1),dummyMap=ALL_context_map,mapLabels=NULL)
p
dev.off()


# plot legend of donor populations
File = paste0('context_plots',directory,'/',version_string,'-DonorLocations-LEGEND.png')
if(withPort) File = gsub('.png','-withPort.png',File)
#png(File,height=8.6*150,width=2.1*150,res=150)
pdf(gsub(".png",".pdf",File),height=(8.6+0.5*withPort+0.5*("donorPopn_318"%in%donorOrder)),width=2.1,bg="transparent")
par(mar=c(0,0,0,0))
plot.new()
legend("center",y.intersp = 1.4,x.intersp = 0.9,title="Group labels",pt.cex=1.4,legend=sapply(donorOrder,getPopnLabel,popnLabels)[1,],bty="n",col="black",pt.bg=recieverLegend2[donorOrder,"colour"],pch=recieverLegend2[donorOrder,"shape"],)
dev.off()

# With sample sizes
#png(File,height=8.6*150,width=2.1*150,res=150)
pdf(gsub(".png","-sampleSizes.pdf",File),height=(8.6+0.5*withPort+0.5*("donorPopn_318"%in%donorOrder)),width=2.5,bg="transparent")
par(mar=c(0,0,0,0))
plot.new()
legendText = sapply(donorOrder,getPopnLabel,popnLabels)[1,]
ssizes = sapply(donorOrder,function(x) sum(DONOR.split.matrix2[,2]==as.numeric(gsub("donorPopn_","",x)) ))
legendText = paste0(legendText," (",ssizes,")")
legend("center",y.intersp = 1.4,x.intersp = 0.9,title="Group labels (size)",pt.cex=1.4,legend=legendText,bty="n",col="black",pt.bg=recieverLegend2[donorOrder,"colour"],pch=recieverLegend2[donorOrder,"shape"],)
dev.off()






###########################
#### Plot coancestry NAfrican vs. sub-Saharan African (with cluster colours). Just raw coancestry.
###########################

# NOTE:
# NAfrican = NorthAfrica.M-A-L      NorthMorocco     WesternSahara     Libya-Algeria           Tunisia    Egypt
# subSaharan = Sub-saharan.Nigeria1 Sub-saharan.Nigeria2         SS.Kenya-LWK         SS.Kenya-MKK             SS.Kenya 

# combine north African coancestry
indsNA = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]%in%theseNA]
indsSS = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]%in%theseSS]
chunkLenthsNA = rowMeans(combinedChunkLengths[rownames(SPAIN.split.matrix),indsNA])
chunkLenthsSS = rowMeans(combinedChunkLengths[rownames(SPAIN.split.matrix),indsSS])

colours = recieverLegend2[paste0("SpainPopn_",SPAIN.split.matrix[names(chunkLenthsNA),2]),"colour"]
shapes = recieverLegend2[paste0("SpainPopn_",SPAIN.split.matrix[names(chunkLenthsNA),2]),"shape"]

# Samples in Portugal-andalucia
classif = (1*SPAIN.split.matrix[,2]==popnLabels$`Spain_Portugal-Andalucia`)

# Compute regression line on all data (assume 0 intercept):
X = chunkLenthsNA[!classif]
Y = chunkLenthsSS[!classif]
regr0 = lm(Y~0+X)  # Intercept at zero
regr = lm(Y~X)  # No intercept at zero <=== More honest!
regr3 = lm(Y~X+I(X^2))  # No intercept at zero <=== More honest!

sregr = summary(regr)
r2 = sregr$r.squared
r20 = summary(regr0)$r.squared

sregr3 = summary(regr3)
r23 = sregr3$r.squared

# get prediction confidence intervals
xx <- data.frame(seq(0, max(chunkLenthsNA), length.out = 30))
colnames(xx) = "X"
pred <- predict.lm(regr, newdata = data.frame(X = xx), interval = "prediction",level=0.99) 

# Compute multiple linear regression line with classification depending whether you are red triangles:
regr2 = lm(chunkLenthsSS~0+chunkLenthsNA*classif)
sregr2 = summary(regr2)
r22 = sregr2$r.squared

######## Scatter plot
#png(paste0('context_plots',directory,'/',version_string,"-meanCoancestryNA-vs-SS.png"),width=1000,height=1000,res=200,bg="transparent")
pdf(paste0('context_plots',directory,'/',version_string,"-meanCoancestryNA-vs-SS.v2.pdf"),width=1000/200,height=1000/200,bg="transparent")

par(cex=0.8,cex.lab=1.2)
plot(chunkLenthsNA,chunkLenthsSS,col="black",pch=shapes,bg=add.alpha(colours,0.8),lwd=0.5,cex=1,
     ylab="mean coancestry with sub-Saharan African donor groups",
     xlab="mean coancestry with north African donor groups")

mtext(side=3,text=bquote(paste(r^2," = ",.(round(r2,3)),sep="")),cex=1)
#mtext(side=3,text=bquote(paste(r^2," = ",.(r2),sep="")),cex=0.7)

#abline(a = coef(regr)[1],b=coef(regr)[2])

#legend("topleft",legend=c(as.expression(bquote(paste(r^2," = ",.(r2),sep=""))),
#                          as.expression(bquote(paste(r^2," = ",.(r22),sep="")))),
#       col=c("red","blue"),cex=1,bty="n",lty=1)

# basic regression
abline(a = coef(regr)[1],b=coef(regr)[2],col="black") #
#mtext(side=4,las=2,at=coef(regr)[1]+coef(regr)[2]*5.4,text=bquote(paste(r^2," = ",.(r2),sep="")),cex=0.7)

#abline(0,1,col="black") # with quadratic term
#myseq =seq(0,8,0.1)
#lines(myseq,coef(sregr3)[1]+coef(sregr3)[2]*myseq+coef(sregr3)[3]*myseq^2)
#mtext(side=4,las=2,at=coef(sregr3)[1]+coef(sregr3)[2]*5.4+coef(sregr3)[3]*5.4^2,text=bquote(paste(r^2," = ",.(r23),sep="")),cex=0.7)

#abline(a = 0,b=coef(regr0)[1],col="black") #
#matlines(xx, pred, lty = c(1, 2, 2), col = "black") 

# two-variable regression
#abline(a = 0 + coef(regr2)[3], b = coef(regr2)[1] + coef(regr2)[4],col="blue") # classif=TRUE
#abline(a = 0 + coef(regr2)[2],b = coef(regr2)[1],col="blue") # classif=FALSE
dev.off()


######## Box plot of ratios

classif2 = SPAIN.split.matrix[names(chunkLenthsNA),2]
classif2factor = factor(classif2,levels=as.numeric(gsub("SpainPopn_","",popOrder)))

rat=(chunkLenthsSS/chunkLenthsNA)  ########  Ratio of SS to NA

#### Check differences in means between group i and everyone else
wil = list()
for( i in as.numeric(levels(classif2factor))){
  wil[[i]] = wilcox.test(formula=rat~(1*classif2==i))
}

#png(paste0('context_plots',directory,'/',version_string,"-meanCoancestryNA-vs-SS-ratioBoxplot.png"),width=1200,height=1000,res=200,bg="transparent")
pdf(paste0('context_plots',directory,'/',version_string,"-meanCoancestryNA-vs-SS-ratioBoxplot.pdf"),width=1200/200,height=1000/200,bg="transparent")
par(opar,mar=c(10,6,2,2),cex=0.6)
boxplot(rat~classif2factor,ylim=range(rat)+c(0,0.01),lwd=0.5,outline=FALSE,names=NA,ylab="Ratio of mean coancestry\n(sub-Subsaharan Africa)/(north Africa) ")
abline(h=mean(rat),lty=3)

for( i in as.numeric(levels(classif2factor))){
    stripchart(rat[classif2==i], at = which(as.numeric(levels(classif2factor))==i), vertical = TRUE,jitter = 0.25, cex=1,xpd=NA,
               add = TRUE, method = "jitter",pch = shapes[classif2==i],axes=FALSE,bg = add.alpha(colours[classif2==i],0.6), col = "black",lwd=0.5)
    pval = wil[[i]]$p.value
    if(pval<0.001) {
      pval=format(pval,scientific = TRUE,digits=2)
      text(x = which(as.numeric(levels(classif2factor))==i),y=max(rat)+0.01-0.003,labels="***")
      } else {
        pval = round(pval,3)
      }
    text(x = which(as.numeric(levels(classif2factor))==i),y=max(rat)+0.01,labels=paste0("p=",pval),xpd=NA,cex=1)
}
# Add labels
nGroups=length(unique(classif2))
#mtext("Iberian clusters",side=1,line=3.5,cex=1.2,col=add.alpha("black",0.5))
xline = strheight("M", units = "user")
par(srt=45)
text(x=1:nGroups+0.1,y=rep(par()$usr[3] - 1.5*xline,nGroups),cex=1.2,pos=2,xpd=NA,labels=sapply(popOrder,getPopnLabel,popnLabels)[1,])
  
dev.off()



