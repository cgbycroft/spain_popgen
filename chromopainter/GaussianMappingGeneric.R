########################################
### Create density maps (based on clustering output)  ####
#######################################
 
DIR="~/Documents/ClareDPhil/DPhil/Spain/densityMaps";#setwd(DIR)
imageHeight=1500
#source("../Other/ClaresFunctions.R")

### get splits and colour information

#plotdir="SPAIN_A2_context"
#serverDIR='SPAIN_A2_context'
#prefix='SPAIN.A2_LOO-SP_PR_HM_NA'
#pre="A2.context"
#chromv=7
#finev='7a'
#version_string = "Spain_A2_context_v7_v7a_contextPainting_v7_DONORpopn.v9" # Spain splits from context fs results
#split.matrix = SPAIN.split.matrix
#resolution="high"
#resolution="low"
#### Which cluster levels do you want to plot? These are columns of the the split.matrix2
#clusterLevels = 1:50 
#clusterLevels = 1:2

###################
#load(paste("../mixtureModel/mixtureModelOutput.",version_string,".Rdata",sep=""),verbose=T)
###################

################## Build spatial grids if needed otherwise just load from file in next step ####################
#source('../densityMaps/getSpatialGrid.R')
##################

################## Get grid information and sample positions

filter = 80
#if( !"distances" %in% ls() ) {
  #if(resolution=="low") load(paste0(DIR,'/Port_A2-SpatialDistances-',filter,'Km-lowRes.Rdata'),verbose=T)
  if(resolution=="low") load(paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-lowRes.Rdata'),verbose=T)
  
  if(resolution=="high"){
    load(paste0(DIR,'/Spain.A2-SpatialDistances-',filter,'Km-highRes.Rdata'),verbose=T)
    distances = distancesZoom
  } 
#}

#### If you want to exclude any samples, do so here...
split.matrix2 = split.matrix
if(sum(grepl("SPAIN_A2",rownames(split.matrix2)))==0) rownames(split.matrix2) = paste0("SPAIN_A2_",rownames(split.matrix2))

#### Exclude samples if necessary
include = !grepl("POPRES",sampleIDs) # always exclude Portugal from Gaussian bit!
#distances = distances[,include]
distances = distances[include,]
sampleIDs = sampleIDs[include] # exclude portuguese from density plot
sampleDistances = sampleDistances[include,include]

#### calculate weights for each cluster level
print("Computing distance from cluster centres...")

splitDist <- sapply(1:ncol(split.matrix2),FUN=function(j){
  s <- unlist(sapply(sampleIDs,FUN=function(x){
               y <- which(sampleIDs==x)
               clust <- split.matrix2[rownames(split.matrix2)==x,j]
               inds <- rownames(split.matrix2)[split.matrix2[,j]==clust]
               indsIndex <- which(sampleIDs%in%inds)
               
               #print(indsIndex)
               if(length(indsIndex)<2) out <- Inf else out <- min(sampleDistances[y,indsIndex[indsIndex!=y]])
               out
              }))
  })


#### Set the minimum distance of a point from it's nearest neighbour in the cluster
print("Evaluating gaussians...")

minDistance=200000

#### Compute Gaussian contributions (takes some minutes)
### Gaussian Parameters ----> this is the final version!!!!!!!!
scaleFactor = 10000
SD=3.5
###
d <- apply(distances/scaleFactor,MARGIN=1,FUN=dnorm,0,SD)
d = t(d) # want rows to be indiivduals and columns to be points in space
sumd <- colSums(d)   # total contribution made to each point (in order to normalise later)
#distancesZoom = distancesZoom/scaleFactor
#dZoom <- apply(distancesZoom,MARGIN=2,dnorm,0,sd=SD)
#sumdZoom <- colSums(dZoom)


#### Set the output directory
dir <- paste(prefix,"_v",chromv,"_v",finev,"-",filter,"Km-clusters-PropTotalGaussianDensity_DIR",sep="")
system(paste('mkdir ',DIR,"/",dir,sep=""))


#### start plotting
print("making density plots...")
if (resolution!="high"){
          SP_grd=sp_grd2
        } 
        if(resolution=="high"){
          SP_grd=sp_grdZoom2
        }
Notinthesea = which(!is.na(over(SP_grd,espMapTotalSansCanProj)[,1]))
dens=d
sumdens <- sumd
      
for (j in clusterLevels){
        
        Weights <- splitDist[,j]    
        Weights[Weights<minDistance] <- 1    
        Weights[Weights>minDistance] <- 1/4  # weight for anything further away than minDistance
        splitSum <- sapply(unique(split.matrix2[,j]),FUN=getSplitSummary,splitColumn=j,split.matrix=split.matrix2,fun="sum",simplify="array",weight=NULL,dens=dens)
        
        colnames(splitSum) <- unique(split.matrix2[,j])
        
        print(j)
        P <- list()
        x=0
        for(i in unique(split.matrix2[,j])){
          splitSumColumn <- which(colnames(splitSum)==as.character(i))
          samples <- rownames(split.matrix2)[which(split.matrix2[,j]==i)]
          
          if(sum(splitSum[,splitSumColumn])>0) {
            
          x = x + 1 # counting how many clusters are being plotted here.
          
          #color <- My.add.alpha(numToColour[2,as.character(i)],alpha=1)
          #color <- recieverLegend2[paste0("SpainPopn_",i),"colour"]
          color <- numToColour[2,as.character(i)]
          #pizza(color)
            #pointCol <- My.add.alpha(color,valueChange=1,alpha=0.1)
            #shape <- numToShape[2,as.character(i%%dim(numToShape)[2])]
            #mappingpoints<- list("sp.points",grand.points.display[which(sampleIDs%in%samples),],col=pointCol,pch=shape,fill=pointCol)
         # sumd2 <- colSums(d[,])
            toPlot <- splitSum[,splitSumColumn]/sumdens
            toPlot[is.na(toPlot)] <- 0
            #toPlot[toPlot<0.02] <- 0
                    
            #p <- plotDenisty(toPlot,color,title="",sp_grd=SP_grd,topAlpha=0.7)
            p <- plotDenisty(toPlot,color,title="",mainMap = "SpainMap_under2Proj",
                             sp_grd=SP_grd,topAlpha=0.7,minDensity = 0,
                             notinthesea = Notinthesea) # use projected map!
            
            if(!is.null(p)) {
              P <- lappend(P,p[[2]])
              print(x)
              if (x==1) Plots <- p[[2]] else Plots <- Plots + p[[2]]
              }
          }
          
        }
    
    if(resolution=="high") assign(paste("highRes",j,sep=""),Plots) else assign(paste("lowRes",j,sep=""),Plots)
    }
    

if(resolution=="high") save(list=ls(pattern="highRes"),file=paste(DIR,"/",dir,"/",gsub("_DIR","",dir),"Gaussian-sd",SD,"-hiResPlots.Rdata",sep="")) else save(list=ls(pattern="lowRes"),file=paste(DIR,"/",dir,"/",gsub("_DIR","",dir),"Gaussian-sd",SD,"-lowResPlots.Rdata",sep=""))   



#### Get positions of samples
stuffForPlots <- getFilteredPointsSpain(maxDist=80000,version="close",zoom=F,split.matrix=split.matrix2,colourVar="Clusters","forPlots",Portugal=TRUE)
grand.points.display <- stuffForPlots$D[rownames(stuffForPlots$D)%in%rownames(split.matrix2),]

### Print plots
print("Printing plots...")
  
for(j in clusterLevels){
  
  File = paste(DIR,"/",dir,"/",gsub("_DIR","",dir),"Gaussian-sd",SD,"-split_",j,"-",resolution,"res-%02d.png",sep="")
  png(File,width=imageHeight*spainAspect,height=imageHeight,res=150)

  Plots <- get(paste(resolution,"Res",j,sep=""))
  #G = update(SpainMap_under2,col.regions="white") +  Plots + update(commapborders2,lwd=0.5,col=hsv(0,0,0.6))
  G = update(SpainMap_under2Proj,col.regions="white") +  Plots + update(commapborders2Proj,lwd=0.5,col=hsv(0,0,0.6)) + update(SpainMap_under2Proj,col.regions="transparent",col="black",lwd=1)

  plotPointsOnMaps2(D=grand.points.display,dend=NULL,plotDendro=F,split.matrix=split.matrix2,pointSize=1,
                    m=j,baseMap=G,minSplitSize=0,title="",
                    dendroPlot=NULL,numToColour=numToColour,numToShapes=numToShape,
                    pointLWD=0.6,zoom=F)
  
  #print(G)
  dev.off()
}

