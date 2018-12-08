########################################
### Run mixture model for SPAIN A2  ####
#######################################

DIR="~/Documents/ClareDPhil/DPhil/Spain/mixtureModel";setwd(DIR)
source("../Other/ClaresFunctions.R")
source("~/Documents/ClareDPhil/DPhil/Spain/mixtureModel/nnls.R")  # Simon's code
source("bootStrapping.R")

load('../chromopainter/SPAIN.A2.badSamples.Rdata')

serverDIR='SP_PR_HM_NA'
prefix='SP_PR_HM_NA'
cont.prefix='PR_HM_NA'
pre="A2"
extra=""
minSpainSplits=10  # unless otherwise stated in version parameters
version_stringStart=""

 {
############# SET VERSION PARAMETERS ############# 
chromv=5 # this is the chromopainter version for the context painting (includes Portugal)
cont.chromv=4  # includes Portugal
SPAIN.splitLevel <- 48
############# 
chromv=6 # this is the chromopainter version for the context painting (excludes Portugal)
cont.chromv=5  # excludes Portugal
SPAIN.splitLevel <- 48
############# 
chromv=7 # this is the chromopainter version for the context painting (includes Portugal)
cont.chromv=7  # includes Portugal
SPAIN.splitLevel <- 48
#############
chromv=8 # this is the chromopainter version for the context painting (excludes Portugal)
cont.chromv=6  # excludes Portugal
SPAIN.splitLevel <- 48
#############
chromv=7 # this is the chromopainter version for the context painting (includes Portugal)
cont.chromv=7  # includes Portugal
SPAIN.splitLevel <- 48
extra = "-noPortDonor"
#############
chromv=7 # this is the chromopainter version for the context painting (includes Portugal)
cont.chromv=7  # includes Portugal
SPAIN.splitLevel <- 48
extra = "-noPortMix"
#############
chromv=7 # this is the chromopainter version for the context painting (includes Portugal)
cont.chromv=7  # includes Portugal
SPAIN.splitLevel <- 48
extra = "-BasqAsSurr"
#############
chromv=7 # this is the chromopainter version for the context painting (includes Portugal)
cont.chromv=7  # includes Portugal
SPAIN.splitLevel <- 48 # try to collapse some Spanish groups 
minSpainSplits=20
#############
chromv=9 # this is the chromopainter version for the context painting (includes Basqu as donor)
cont.chromv=8  # includes Portugal and includes Basques as donor
SPAIN.splitLevel <- 48 # try to collapse some Spanish groups 
minSpainSplits=20
#############
chromv=10 # this is the chromopainter version for the context painting (more Irish)
cont.chromv=9  # includes Portugal (more Irish)
SPAIN.splitLevel <- 48 # try to collapse some Spanish groups but only in Galicia
minSpainSplits=20
extra="GalCombine"
donorVersion=9 
#############
chromv=10 # this is the chromopainter version for the context painting (more Irish)
cont.chromv=9  # includes Portugal (more Irish)
SPAIN.splitLevel <- 48 # try to collapse some Spanish groups but only in Galicia
minSpainSplits=20
extra="-noPortDonor"
donorVersion=9 
#############
chromv=10 # this is the chromopainter version for the context painting (more Irish)
cont.chromv=9  # includes Portugal (more Irish)
SPAIN.splitLevel <- 48 # try to collapse some Spanish groups but only in Galicia
minSpainSplits=20
extra="-noPortMix"
donorVersion=9 
#############
chromv=6 # this is the chromopainter version for the context painting (more Irish)
cont.chromv=6  # includes Portugal (more Irish)
cont.prefix='LOO-SP_PR_HM_NA'
SPAIN.splitLevel <- 48 # try to collapse some Spanish groups but only in Galicia
minSpainSplits=20
extra=""
donorVersion=9
version_stringStart= "v1_v7a_48_minSize20_contextPainting_v10_v3a-DONORpopn.v9-sumsGalCombine"
version_string='Spain_v1_v7a_contextPainting_v6_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
#############
cont.prefix='LOO-SP_PR_HM_NA'  #exclude Italy,France and Portugal from everything
chromv=6
version_stringStart='Spain_v1_v7a_contextPainting_v6_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
extra="-noItalyFrancePortugal"
version_string=paste(version_stringStart,extra,sep="")
#############
cont.prefix='LOO-SP_PR_HM_NA'  # exclude portugal from everything
chromv=6
version_stringStart='Spain_v1_v7a_contextPainting_v6_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
extra="-noPortMix"
version_string=paste(version_stringStart,extra,sep="")
}

#############
chromv=7 # this is the chromopainter version for the context painting with portuguese as recipients
cont.chromv=7  # includes Portugal as recipient
cont.prefix='LOO-SP_PR_HM_NA'
extra="portAsRecipient"
donorVersion=9
version_stringStart= 'Spain_v1_v7a_contextPainting_v6_DONORpopn.v9'
version_string='Spain_v1_v7a_contextPainting_v7_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
#############
chromv=6 # this is the chromopainter version for the context painting (more Irish)
cont.chromv=6  # includes Portugal (more Irish)
cont.prefix='LOO-SP_PR_HM_NA'
SPAIN.splitLevel <- 14 # pick a high part in the tree and set min size to 20
minSpainSplits=20
extra="-Spain14Clusters"
donorVersion=9
version_string=paste('Spain_v1_v7a_contextPainting_v6_DONORpopn.v9',extra,sep="")
#############
chromv=7 # this is the chromopainter version for the context painting with portuguese as recipients
cont.chromv=7  # includes Portugal as recipient
cont.prefix='LOO-SP_PR_HM_NA'
extra="Spain14Clusters-portAsRecipient"
SPAIN.splitLevel <- 14 # pick a high part in the tree and set min size to 20
minSpainSplits=20
donorVersion=9
version_string=paste('Spain_v1_v7a_contextPainting_v7_DONORpopn.v9',extra,sep="")
#############
chromv=7 # this is the chromopainter version for the context painting with portuguese as recipients
cont.chromv=7  # includes Portugal as recipient
cont.prefix='LOO-SP_PR_HM_NA'
extra="-noItaly"  # exclude italy entirely (ie. re-normalise copying matrix)
minSpainSplits=20
donorVersion=9
version_stringStart='Spain_v1_v7a_contextPainting_v7_DONORpopn.v9'
version_string=paste(version_stringStart,extra,sep="")
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
#############
chromv=7 # this is the chromopainter version for the context painting with portuguese as recipients
cont.chromv=7  # includes Portugal as recipient
cont.prefix='LOO-SP_PR_HM_NA'
extra="-noItalyDonor"  # exclude italy as surrogate in the mixture for Spain (but keep as an item in copy vector)
minSpainSplits=20
donorVersion=9
version_stringStart='Spain_v1_v7a_contextPainting_v7_DONORpopn.v9'
version_string=paste(version_stringStart,extra,sep="")
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
#############
chromv=9 # this is the chromopainter version for the context painting with portuguese as recipients
cont.chromv=9  # includes Portugal as recipient + basques as donors
cont.prefix='LOO-SP_PR_HM_NA'
extra=""
minSpainSplits=20
donorVersion=10
version_stringStart='Spain_v1_v7a_contextPainting_v7_DONORpopn.v9'
version_string='Spain_v1_v7a_contextPainting_v9_DONORpopn.v10'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""))
#############
chromv=9 # this is the chromopainter version for the context painting with portuguese as recipients
cont.chromv=9  # includes Portugal as recipient + basques as donors
cont.prefix='LOO-SP_PR_HM_NA'
extra="-noItalyDonor"  # exclude italy as surrogate in the mixture for Spain (but keep as an item in copy vector)
minSpainSplits=20
donorVersion=10
version_stringStart='Spain_v1_v7a_contextPainting_v9_DONORpopn.v10'
version_string=paste(version_stringStart,extra,sep="")
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=TRUE)
#############
# this is the version with donorVersion 9, but Spain target popns based on spain context fs run v7_v7a.
# NOTE Jan 2016: I haven't actually run the mixture model for this yet, just updated the labels
version_stringStart='Spain_v1_v7a_contextPainting_v7_DONORpopn.v9'
donorVersion=9
cont.prefix='LOO-SP_PR_HM_NA'
cont.chromv=7
version_string='Spain_A2_context_v7_v7a_contextPainting_v7_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=T)
SPAIN.split.matrix = read.table('../chromopainter/SPAIN_A2_context/SPAIN.split.matrix_A2.context.v7_v7a.txt',header=T,stringsAsFactors=F)
SPAIN.splitLevel = "contextClust"
#############
# as above, but version 10a
version_stringStart='Spain_A2_context_v7_v7a_contextPainting_v7_DONORpopn.v9'
donorVersion=9
cont.prefix='LOO-SP_PR_HM_NA'
cont.chromv=7
version_string='Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=T)
recieverLegendOrig = recieverLegend2
SPAIN.split.matrix = read.table('../chromopainter/SPAIN_A2_context/SPAIN.split.matrix_A2.context.v7_v10a.txt',header=T,stringsAsFactors=F)
SPAIN.splitLevel = "contextClust"
load(paste0("../chromopainter/SPAIN.A2_LOO-SP_PR_HM_NA.ProcessedResults_v1-v7_v10a.RData"),verbose=TRUE)
#############
# as above but exclude high france outliers
version_stringStart='Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=T)
extra = "-excHighBasFrance"
cont.chromv=7
cont.prefix='LOO-SP_PR_HM_NA'
version_string=paste0(version_stringStart,extra)
#load(paste("../mixtureModel/mixtureModelOutput.",version_string,".Rdata",sep=""),verbose=T)
#############
# as above but use means instead of sums!
version_string = "Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9-excHighBasFrance-means"
version_stringStart = "Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9-excHighBasFrance"
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=T)
model="means"
#############
# as three above, but with basques as donors
version_stringStart='Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=T)
SPAIN.split.matrix1 = SPAIN.split.matrix
donorVersion=12
load(paste('../mixtureModel/donorPopns/DONOR.popn.objects_v',donorVersion,'.Rdata',sep=""),verbose=TRUE)
cont.prefix='LOO-SP_PR_HM_NA'
cont.chromv=10
load(paste0("../chromopainter/SPAIN.A2.ProcessedResults_v1-Spain_v1_v7a.new.RData"),verbose=TRUE)
SPAIN.splitLevel = "Split_27"
version_string='Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12'
#############
# as three above, but Combining Galicia Coast into one group as in GT analysis
version_stringStart='Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9'
load(paste("../mixtureModel/mixtureModelOutput.",version_stringStart,".Rdata",sep=""),verbose=T)
SPAIN.split.matrix1 = SPAIN.split.matrix
donorVersion=12
load(paste('../mixtureModel/donorPopns/DONOR.popn.objects_v',donorVersion,'.Rdata',sep=""),verbose=TRUE)
cont.prefix='LOO-SP_PR_HM_NA'
cont.chromv=10
load(paste0("../chromopainter/SPAIN.A2.ProcessedResults_v1-Spain_v1_v7a.new.RData"),verbose=TRUE)
SPAIN.splitLevel = "Split_27"
version_string='Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine'
#############

#model="sums"

#donorVersion=1   # note: donors must at least contain all the donors in the chromv version - check donorPopn.R
#donorVersion=3
#donorVersion=4  # includes more popres samples
#donorVersion=5  # combines swiss groups
#donorVersion=6  # as 5 but with Basque included as donor
#donorVersion=8  # as 5 but with differnet NA splits
#donorVersion=9  # as 8 but with more Irish and differnet NA splits

spain.chromv=1
spain.finev='7a.new'
cont.finev='3a'

SPAIN.splitLevel2 <- SPAIN.splitLevel

#version_string <- paste('v',spain.chromv,'_v',spain.finev,'_',SPAIN.splitLevel,'_contextPainting_v',chromv,'_v',cont.finev,'-DONORpopn.v',donorVersion,'-',model,extra,sep="")
#if((version_stringStart!="")) {
#  version_string <- paste('v',spain.chromv,'_v',spain.finev,'_',SPAIN.splitLevel,'_minSize',minSpainSplits,'_contextPainting_v',chromv,'_v',cont.finev,'-DONORpopn.v',donorVersion,'-',model,extra,sep="")
#} else {
#}
################# get some donorPopn Info and colour things for legends

colourSet <- c("yellow","red","blue","black")

#load(paste('../samples/numToColourForGaussian_A2_v',spain.chromv,'_v',spain.finev,'.Rdata',sep=""))
#load(paste('../mixtureModel/donorPopns/DONOR.popn.objects_v',donorVersion,'.Rdata',sep=""),verbose=TRUE)

############### INPUT chromopainter/finestructure data (if needed)###############
  
if(cont.prefix=='LOO-SP_PR_HM_NA'){
  load(paste('../chromopainter/',cont.prefix,".ChromopainterOutput_v",cont.chromv,".Rdata",sep=""))  
  load(paste('../chromopainter/new.split.matrix.',pre,'_v',spain.chromv,'_v',spain.finev,'.Rdata',sep=""))
  load(paste('../chromopainter/SPAIN.',pre,'.TreeList_v',spain.chromv,'_v',spain.finev,'.Rdata',sep=""))

    SPAIN.split.matrix_raw <- get(paste('new.split.matrix.',pre,'_v',spain.chromv,'_v',spain.finev,sep=""))
  
} else {
    #a <- paste("../chromopainter/chromo_out/",prefix,"_v",chromv,".renamed.chunklengths.out",sep="")
    a <- grep(paste(prefix,"_v",chromv,".*\\.chunklengths.out",sep=""),list.files("../chromopainter/chromo_out/"),value=T,perl=T)
    SPAINchunklengths_raw=fixNames(as.matrix(read.table(paste("../chromopainter/chromo_out/",a,sep=""),row.names=1,header=T)))  ## chromopainter chunkcounts file for context painting
      rownames(SPAINchunklengths_raw)[which(rownames(SPAINchunklengths_raw)=="SPAIN_A2_150X")] <- "SPAIN_A2_150XX"
    
    #a <- paste("../chromopainter/chromo_out/",prefix,"_v",chromv,".renamed.chunkcounts.out",sep="")
    a <- grep(paste(prefix,"_v",chromv,".*\\.chunkcounts.out",sep=""),list.files("../chromopainter/chromo_out/"),value=T,perl=T)
    SPAINchunkcounts_raw=fixNames(as.matrix(read.table(paste("../chromopainter/chromo_out/",a,sep=""),row.names=1,header=T)))  ## chromopainter chunkcounts file for context painting
      rownames(SPAINchunkcounts_raw)[which(rownames(SPAINchunkcounts_raw)=="SPAIN_A2_150X")] <- "SPAIN_A2_150XX"
    
    SPAINmeanchunklengths_raw=SPAINchunklengths_raw/SPAINchunkcounts_raw
    
    
    #b <- paste("../chromopainter/chromo_out/",cont.prefix,"_v",cont.chromv,".renamed.chunklengths.out",sep="")
    b <- paste("../chromopainter/chromo_out/",cont.prefix,"_v",cont.chromv,".chunklengths.out",sep="")
    DONORchunklengths_raw=fixNames(as.matrix(read.table(b,row.names=1,header=T)))   ## chromopainter chunklengths file
    
    #b <- paste("../chromopainter/chromo_out/",cont.prefix,"_v",cont.chromv,".renamed.chunkcounts.out",sep="")
    b <- paste("../chromopainter/chromo_out/",cont.prefix,"_v",cont.chromv,".chunkcounts.out",sep="")
    DONORchunkcounts_raw=fixNames(as.matrix(read.table(b,row.names=1,header=T)))   ## chromopainter chunkcounts file
    
    DONORmeanchunklengths_raw=DONORchunklengths_raw/DONORchunkcounts_raw
    
    #checks on split matrices

    if (!setequal(colnames(DONORchunklengths_raw),colnames(SPAINchunklengths_raw))) {
      stop("ERROR: Spain and Context chromopainter versions do not correspond. check version parameters")
      print('hello')
    }

    load(paste('../chromopainter/new.split.matrix.',pre,'_v',spain.chromv,'_v',spain.finev,'.Rdata',sep=""))
    load(paste('../chromopainter/SPAIN.',pre,'.TreeList_v',spain.chromv,'_v',spain.finev,'.Rdata',sep=""))

    SPAIN.split.matrix_raw <- get(paste('new.split.matrix.',pre,'_v',spain.chromv,'_v',spain.finev,sep=""))

}


# read in data for bootsrapping

#SCP(filename=paste('chromopainter_files/',serverDIR,'/v',chromv,'/combined2/byChrom/*.chunklengths.out',sep=""),"from",recievedir=paste('mixtureModel/bootstrapData',sep=""),spain=T)
SCP(filename=paste('chromopainter_files/',cont.prefix,'/v',cont.chromv,'/combined2/byChrom/*.chunklengths.out',sep=""),"from",recievedir=paste('mixtureModel/bootstrapData',sep=""),spain=T)

if(version_stringStart==""){
bad <- badSamples[spain.finev]
    SPAIN.split.matrix <- excludeSplits(SPAIN.split.matrix_raw,exclusionSplits=bad)
    exclusion_samples <- rownames(SPAIN.split.matrix_raw)[which(SPAIN.split.matrix_raw[,bad]==bad)]  
    list <- get(paste('TreeList_v',spain.chromv,'_v',spain.finev,sep=""))
    Tree =  drop.tip(list$ttree,tip=exclusion_samples)
    SpainTree = Tree
#    dend <- dendro_data(myapetodend(Tree,factor=1))
}

if(donorVersion==10){
bad <- badSamples[spain.finev]
basqSamples=rownames(SPAIN.split.matrix)[SPAIN.split.matrix[,2]%in%c(3,15,29)]
    SPAIN.split.matrix <- SPAIN.split.matrix[!rownames(SPAIN.split.matrix)%in%basqSamples,]
    Tree =  drop.tip(SpainTree,tip=basqSamples)
    SpainTree = Tree
#    dend <- dendro_data(myapetodend(Tree,factor=1))
}

if(donorVersion==12){
    basqSamples=paste0("SPAIN_A2_",rownames(SPAIN.split.matrix)[SPAIN.split.matrix[,2]%in%c(3,18)])
    sm <- SPAIN.split.matrix1[!rownames(SPAIN.split.matrix1)%in%basqSamples,]
    SpainTree$tip.label = paste0("SPAIN_A2_",SpainTree$tip.label)
    SpainTree =  drop.tip(SpainTree,tip=basqSamples)
    # replace split.matrix values in matrix
    rownames(SPAIN.split.matrix) =paste0("SPAIN_A2_",rownames(SPAIN.split.matrix))
    notBasq = rownames(SPAIN.split.matrix)[!rownames(SPAIN.split.matrix)%in%basqSamples]
    sm[notBasq,] = SPAIN.split.matrix[notBasq,]
    port = grep("POPRES",rownames(sm),value=TRUE)
    sm[port,1] = NA; sm[port,2] = 100
    colnames(sm) = colnames(SPAIN.split.matrix)
    SPAIN.split.matrix = sm
}


################# Summarize and combine matrices ######################

# combine chunklength matrices

if(cont.prefix == 'LOO-SP_PR_HM_NA'){
    combinedChunkLengths <- get(paste("chunklengthsmatrix_v",cont.chromv,sep=""))
    combinedMeanChunkLengths <- get(paste("meanchunklengthsmatrix_v",cont.chromv,sep=""))
    combinedChunkLengths[is.na(combinedChunkLengths)]=0
    combinedMeanChunkLengths[is.na(combinedMeanChunkLengths)]=0
    
}else{  
    combinedChunkLengths <- rbind(SPAINchunklengths_raw,DONORchunklengths_raw)
    combinedMeanChunkLengths <- rbind(SPAINmeanchunklengths_raw,DONORmeanchunklengths_raw)  
}

# set split level for SPAIN

if(version_stringStart==""){

SPAIN.splitLevel <- paste("Split",SPAIN.splitLevel,sep="_")
SPAIN.splits <- unique(SPAIN.split.matrix[,SPAIN.splitLevel])

# combine galician groups only
if(extra%in%c("GalCombine","-Spain14Clusters","Spain14Clusters-portAsRecipient")){
  Gal <- which(SPAIN.split.matrix[,2]==2)
  SPAIN.split.matrix[,1] <- SPAIN.split.matrix[,SPAIN.splitLevel]
  SPAIN.split.matrix[,2] <- SPAIN.split.matrix[,1]
  SPAIN.split.matrix[Gal,2] <- SPAIN.split.matrix[Gal,14]
  SPAIN.split.matrix <- SPAIN.split.matrix[,1:2]
  colnames(SPAIN.split.matrix) <- c(paste(SPAIN.splitLevel,"_original",sep=""),SPAIN.splitLevel)
  
} else {

# collapse small groups in SPAIN
  x <- SPAIN.split.matrix[,SPAIN.splitLevel]
  names(x) <- rownames(SPAIN.split.matrix)
  newSPAINsplits <- collapseSplits(x,minSplitSize=minSpainSplits,Tree=Tree)
  SPAIN.splits <- unique(newSPAINsplits)
  SPAIN.split.matrix <- SPAIN.split.matrix[,1:2]
  SPAIN.split.matrix[,1] <- x
  SPAIN.split.matrix[,2] <- newSPAINsplits
  colnames(SPAIN.split.matrix) <- c(paste(SPAIN.splitLevel,"_original",sep=""),SPAIN.splitLevel)

}

# set split level for CONTEXT
splitLevel <- "new"
DONOR.splitLevel <- paste("Split",splitLevel,sep="_")
Port = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==6]
}

##### Extra things ######

# include Basque as surrogate, combine Basque1 populations: 29, 15, 3 into one
if((spain.chromv==1)&(spain.finev=="7a")){
  Basq <- c(29,15,3)
  BasqNum <- as.numeric(paste(Basq,collapse=""))
}

if(extra=="-BasqAsSurr"){  
  SPAIN.splits <- c(SPAIN.splits[!SPAIN.splits%in%Basq],BasqNum)
  SPAIN.split.matrix[SPAIN.split.matrix[,2]%in%Basq,2] <- BasqNum
}

if(version_string=='Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine'){
  # combine Galicia coast into one group (call it number 2). This is what I did for the GT analysis v16.
  
  SPAIN.split.matrix[SPAIN.split.matrix[,2]%in%unlist(SPAINpopnLabels[names(SPAINpopnLabels)=="Spain_Galicia_coast"]),]=2
  
}


# donor populations
# if non-Portugal version then remove samples in that donor population from mixture
if(chromv==8){
port <- rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==6]
combinedChunkLengths <- combinedChunkLengths[!rownames(combinedChunkLengths)%in%port,!colnames(combinedChunkLengths)%in%port]
combinedMeanChunkLengths <- combinedMeanChunkLengths[!rownames(combinedMeanChunkLengths)%in%port,!colnames(combinedMeanChunkLengths)%in%port]
}

# exclude portuguese samples from all mixture
if((extra=="-noPortMix")){
  port <- which(DONOR.split.matrix[,2]==6)
  DONOR.split.matrix <- DONOR.split.matrix[-port,]
}

if(extra=="-noItalyFrancePortugal"){
  extraPopns=unlist(popnLabels[grep("France|Italy|Portugal",names(popnLabels))])
  toRemove <- which(DONOR.split.matrix[,2]%in%extraPopns)
  DONOR.split.matrix <- DONOR.split.matrix[-toRemove,]  
}

if(extra=="-noItaly"){
  extraPopns=unlist(popnLabels[grep("Italy",names(popnLabels))])
  toRemove <- which(DONOR.split.matrix[,2]%in%extraPopns[1])
  DONOR.split.matrix <- DONOR.split.matrix[-toRemove,]  
}
  
if(grepl("portAsRecipient",extra)){
  Port=which(DONOR.split.matrix[,2]==6)
  newPort=DONOR.split.matrix[Port,]
  colnames(newPort )=colnames(SPAIN.split.matrix)
  newPort[,2]=max(SPAIN.split.matrix[,2])+1
  newPort[,2]=49
  SPAIN.split.matrix=rbind(SPAIN.split.matrix,newPort)
  DONOR.split.matrix=DONOR.split.matrix[-Port,]
}

if(extra=="-excHighBasFrance"){
  fra = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==17]
  bas = rownames(SPAIN.split.matrix)[SPAIN.split.matrix[,2]%in%c(popnLabels$Spain_Basque1,popnLabels$Spain_Basque2)]
  toTest = combinedChunkLengths[,fra]
  hist(colMeans(as.matrix(toTest[bas,])),breaks=50)
  HighBasFrance = colnames(toTest)[colMeans(as.matrix(toTest[bas,])) > 8]
 # exclude from DONOR.split.matrix
  DONOR.split.matrix =  DONOR.split.matrix[!rownames(DONOR.split.matrix)%in%HighBasFrance,]
  # checking after the fact 
  bas2 = rownames(S.split.matrix)[S.split.matrix[,2]%in%c(SPAINpopnLabels$Spain_Basque1,SPAINpopnLabels$Spain_Basque2)]
  bas2 = paste0("SPAIN_A2_",bas2)
  toTest2 = combinedChunkLengths[bas2,fra]
  hist(colMeans(as.matrix(toTest2)),breaks=50)
  hist(apply(as.matrix(toTest2),2,max),breaks=50)
  quantile(colMeans(as.matrix(toTest2)))
}


DONOR.split.matrix <- DONOR.split.matrix[which(rownames(DONOR.split.matrix)%in%colnames(combinedChunkLengths)),]
DONOR.splits <- unique(DONOR.split.matrix[,DONOR.splitLevel])

  
############### create summary matrix

DONOR.splits <- unique(DONOR.split.matrix[,2])
SPAIN.splits <- unique(SPAIN.split.matrix[,2])

donorPopns <- paste("donorPopn_",DONOR.splits,sep="")
SpainPopns <- paste("SpainPopn_",SPAIN.splits,sep="")

ALLpopns <- c(SpainPopns,donorPopns)
ALLinds <- rownames(combinedChunkLengths)

combinedChunkLengths = as.data.frame(combinedChunkLengths)
combinedMeanChunkLengths = as.data.frame(combinedMeanChunkLengths)

summaryOutChunkLengths <- getCountSummaries(combinedChunkLengths,DONOR.splits,SPAIN.splits,
         DONOR.split.matrix,SPAIN.split.matrix,
         DONOR.splitLevel,SPAIN.splitLevel,
         ALLpopns,ALLinds)

SUMMARYmat_1 <- summaryOutChunkLengths[[1]]
SUMMARYmat_1mean <- summaryOutChunkLengths[[2]]
SUMMARYmat <- summaryOutChunkLengths[[3]]
SUMMARYmat_mean <- summaryOutChunkLengths[[4]]

summaryOutMeanChunkLengths <- getCountSummaries(combinedMeanChunkLengths,DONOR.splits,SPAIN.splits,
         DONOR.split.matrix,SPAIN.split.matrix,
         DONOR.splitLevel,SPAIN.splitLevel,
         ALLpopns,ALLinds)

# for doing donor groups!

# if version without portugal as donor, but allow portugal copying then remove portugal from donor list
if(extra=="-noPortDonor"){
  port <- "donorPopn_6"
  donorPopns <- donorPopns[-which(donorPopns==port)]
}

# if version with Basque as a surrogate in model
if(extra=="-BasqAsSurr"){
    donorPopns <- c(donorPopns,paste("SpainPopn",BasqNum,sep="_"))
}

if(extra=="-noItalyDonor"){
  It <- "donorPopn_5"
  donorPopns <- donorPopns[-which(donorPopns==It)]
}

#png('test.png',width=1000,height=1000,res=100)
#toPlot <- SUMMARYmat
#    plotFinestructure(toPlot,labelsx=rownames(toPlot),labelsy=colnames(toPlot))
#dev.off()

#system(paste('mkdir ',spainDIR,'/mixtureModel/context_plots/',filename,'_DIR',sep=""))
#filename=paste('SPAIN.A2.Context.pre-mix_',version_string,sep="")
#png(paste('context_plots/',filename,'_DIR/',filename,'-%02d.png',sep=""),width=500,res=100)
#    for (i in c(1:length(SpainPopns))){
#    par(las=2)
#      b <- barplot(SUMMARYmat[i,]/sum(SUMMARYmat[i,]),main=rownames(SUMMARYmat)[i],ylim=c(0,1))
    #print(b)
#    }
#dev.off()


################# Run mixture model #####################

# NOTE: here 'donorPopns' are those that are included as covariates in the mixture (not necessarily the same as the list of popns in the copying vectors)
newmat <- SUMMARYmat[donorPopns,]
mixmat <- matrix(0,nrow=nrow(SUMMARYmat),ncol=nrow(newmat))
  colnames(mixmat) <- donorPopns
  rownames(mixmat) <- rownames(SUMMARYmat)
mixmat_output <- list()

for(pop in ALLpopns){
  print(pop)
  tmpmat <- newmat[rownames(newmat)!=pop,colnames(newmat)!=pop] # the rows of this are the copying vectors which are the input vectors X for the nnls. 

  # treat donor popns separately (i.e need to recompute SUMMARYmat without copying from its own population)
  if( !grepl("donor",pop) ){
    ourmixALL <- getoverallfit(tmpmat/rowSums(tmpmat),
                        SUMMARYmat[pop,colnames(newmat)!=pop]/sum(SUMMARYmat[pop,colnames(newmat)!=pop]))
     
  } else {
    predMat = donorSummaries[[pop]]$predMat
    toFit = donorSummaries[[pop]]$toFit
    ourmixALL <- getoverallfit(predMat,toFit)    
  }
  
  mixmat_output[[pop]] <- ourmixALL
  ourmix <- ourmixALL$x
  ourmix <- ourmix[ourmix>0]
 # ourmix <- ourmix/sum(ourmix)
  mixmat[pop,colnames(mixmat)%in%names(ourmix)] <- ourmix
}

newmat <- SUMMARYmat[donorPopns,]
#SUMMARYmat_inds <- rbind(SUMMARYmat_1[grep("SPAIN_A2",rownames(SUMMARYmat_1)),],newmat)
SUMMARYmat_inds <- rbind(SUMMARYmat_1[rownames(SUMMARYmat_1)%in%rownames(SPAIN.split.matrix),],newmat)

mixmat_inds <- matrix(0,nrow=nrow(SUMMARYmat_inds),ncol=nrow(newmat))
  colnames(mixmat_inds) <- donorPopns
  rownames(mixmat_inds) <- rownames(SUMMARYmat_inds)
mixmat_output_inds <- list()

for(pop in rownames(mixmat_inds)){
  print(pop)
  tmpmat <- newmat[rownames(newmat)!=pop,colnames(newmat)!=pop]
  ourmixALL <- getoverallfit(tmpmat/rowSums(tmpmat),
                        SUMMARYmat_inds[pop,colnames(newmat)!=pop]/sum(SUMMARYmat_inds[pop,colnames(newmat)!=pop]))
  mixmat_output_inds[[pop]] <- ourmixALL
  ourmix <- ourmixALL$x
  ourmix <- ourmix[ourmix>0]
 # ourmix <- ourmix/sum(ourmix)
  mixmat_inds[pop,colnames(mixmat_inds)%in%names(ourmix)] <- ourmix
}


## Bootstrap - sample with replacement from each cluster and recalculate mixmat
# input SPAINnorm and DONORnorm

#### new version

# read in bootstrap raw data

print('reading in bootstrap files')
inputDir=paste(DIR,'/bootstrapData',sep="")
summaryType ='chunklengths'
seq=c(1:22)

file = paste(version_string,'-bootstrapData.Rdata',sep="")
if (!file%in%list.files(paste(DIR,'/bootstrapData',sep=""))){

  
for (i in seq){
    print(i)
    #File <- list.files(inputDir,pattern=paste(cont.prefix,'_v',cont.chromv,'.chr',i,'.',summaryType,'.out',sep=""))
    File <- list.files(inputDir,pattern=paste(cont.prefix,'.chr',i,'_v',cont.chromv,'_Popn_SPAIN_A2.',summaryType,'.out',sep=""))    
    print(File)
    data <- read.table(paste(inputDir,'/',File,sep=""),header=T,stringsAsFactors=F)
    
    if(i==1) {
      RecipientOrder = data[,1]
      DonorOrder <- colnames(data)[-1]
      Inds <- array(dim=c(length(seq),nrow(data),ncol(data)-1))
    }
    
    thisDataOrderRows = data[,1]
    thisChrom = as.matrix(data[,-1])
    
    print( sum(thisDataOrderRows!=RecipientOrder) )
    print( sum(colnames(thisChrom)!=DonorOrder) )
    Inds[i,,] <- thisChrom[match(RecipientOrder,thisDataOrderRows),DonorOrder]
}

#test = apply(Inds,3,colSums)
      
#RecipientOrder <- data[,1]
#DonorOrder <- colnames(data)[-1]

save(Inds,RecipientOrder,DonorOrder,file=paste(DIR,'/bootstrapData/',version_string,'-bootstrapData.Rdata',sep=""))

} else {
  load(paste(DIR,'/bootstrapData/',version_string,'-bootstrapData.Rdata',sep=""),verbose=TRUE)
}
  
#Spain.Popns <- SPAIN.split.matrix[match(RecipientOrder,paste("SPAIN_A2",rownames(SPAIN.split.matrix),sep="_")),2,drop=F]
Spain.Popns <- SPAIN.split.matrix[match(RecipientOrder,rownames(SPAIN.split.matrix)),2,drop=F]

# get bootstrap samples (dim=Nrecipients x Ndonors)

File = paste(version_string,'-bootstrapDataPostBoot.Rdata',sep="")
if(!File%in%list.files(paste(DIR,'/bootstrapData',sep=""))){
  bootOut <- lapply(c(1:1000),FUN=getOneBootstrap,Inds=Inds,popnList=Spain.Popns)
  #save(bootOut,file=paste(DIR,'/bootstrapData/',File,sep=""))   <=== files are 10s GB! 
} else {
  #load(paste(DIR,'/bootstrapData/',File,sep=""))
}
   
# summarise each by popns and run nnls
bootstrappedmats <- lapply(bootOut,FUN=runNNLS,SUMmat=SUMMARYmat_1,Y=SUMMARYmat)

# new version to avoid large dataset!
 
############
# old version
bootstrappedmats <- lapply(c(1:100),FUN=mixbootstrap,SUMmat=SUMMARYmat_1,Y=SUMMARYmat)
############

meanbootstrap <- bootstrappedmats[[1]]
percentile95bootstrap <- bootstrappedmats[[1]]
percentile05bootstrap <- bootstrappedmats[[1]]

for (i in c(1:nrow(bootstrappedmats[[1]]))) {
        for (j in c(1:ncol(bootstrappedmats[[1]]))) {
             meanbootstrap[i,j] <- mean(sapply(bootstrappedmats,FUN=function(x) x[i,j]))
             #percentile95bootstrap[i,j] <- quantile(sapply(bootstrappedmats,FUN=function(x) x[i,j]),0.95)
             percentile95bootstrap[i,j] <- quantile(sapply(bootstrappedmats,FUN=function(x) x[i,j]),0.975)
             #percentile05bootstrap[i,j] <- quantile(sapply(bootstrappedmats,FUN=function(x) x[i,j]),0.05)
             percentile05bootstrap[i,j] <- quantile(sapply(bootstrappedmats,FUN=function(x) x[i,j]),0.025)

        }
}


################################## POPN labels (manual)

if(version_stringStart==""){
    
  if((spain.chromv==1)&(spain.finev=="7a")){
      SPAINpopnLabels <- list(c(2,37,19,38,5,28,25,17,32,26,23,36,11,41,14,39,9,27,30,4,24,22,10,42,46,44,49),c(12,31),c(20),c(8,35),c(3,29,15,29153),c(7,47),c(13,33),c(34),c(40),c(21),c(43),c(6),c(1),c(18),c(48))
      names(SPAINpopnLabels) <- c("Spain_Galicia_coast","Spain_Galicia_inland","Spain_Asturias","Spain_Cataluna","Spain_Basque1","Spain_Basque2","Spain_Aragon-Valencia","Spain_LaRioja","Spain_Murcia","Spain_SouthCoast","Spain_Valencia","Spain_central","Spain_west","Spain_northCentral","Spain_Baleares")
  }
  popnLabels <- append(DONORpopnLabels,SPAINpopnLabels)
  
  ### augment recieverLegend with Spain samples and other labels ###
  BasqNum=29153
  numToColour <- cbind(numToColour,c(BasqNum,numToColour[2,"3"]))
                       colnames(numToColour)[numToColour[1,]==as.character(BasqNum)] <- as.character(BasqNum)
  AllSpainNos <- unlist(SPAINpopnLabels)
  spainLegend <- cbind(paste("SpainPopn",AllSpainNos,sep="_"),as.vector(numToColour[2,as.character(AllSpainNos)]),
                       numToShape[2,as.character(AllSpainNos%%dim(numToShape)[2])])
      rownames(spainLegend) <- spainLegend[,1]
      colnames(spainLegend) <- colnames(DONORLegend)
  
  recieverLegend <- rbind(DONORLegend,spainLegend)
  
  popnLabelsCut <- sapply(popnLabels , "[[", 1)
  newDon <- recieverLegend[paste("donorPopn",popnLabelsCut[grep("Spain",names(popnLabels),invert=T)],sep="_"),]
  newSp <- recieverLegend[paste("SpainPopn",popnLabelsCut[grep("Spain",names(popnLabels))],sep="_"),]
  new <- rbind(newDon,newSp)
      rownames(new) <- new$Factor <- names(popnLabels)
  recieverLegend2 <- rbind(recieverLegend,new)
  recieverLegend2$shape <- as.numeric(recieverLegend2$shape)

}

if(grepl("portAsRecipient",extra)){
  popnLabels$Spain_Portugal = 49
  popnLabels$Spain_Galicia_coast = popnLabels$Spain_Galicia_coast[!popnLabels$Spain_Galicia_coast==49]
}

if(version_string=="Spain_A2_context_v7_v7a_contextPainting_v7_DONORpopn.v9"){
  popnLabels = popnLabels[!grepl("Spain",names(popnLabels))]
  popnLabels[["Spain_south-Portugal"]] = c(1)
  popnLabels[["Spain_High_SubSaharan"]] = c(2)
  popnLabels[["Spain_central"]] = c(3)
  popnLabels[["Spain_Cataluna-Aragon"]] = c(4)
  popnLabels[["Spain_Basque"]] = c(5)
  popnLabels[["Spain_Galicia-Portugal"]] = c(6)
  popnLabels[["Spain_west"]] = c(7)
  
  splits = unique(SPAIN.split.matrix[,2])
  
  newLegend = makeFactorLegend(splits,shapes=numToShape[2,as.character(splits%%5)],newColours=numToColour[2,splits])
  newLegend["2","colour"] = "purple"
  rownames(newLegend) = newLegend$Factor = paste0("SpainPopn_",rownames(newLegend)) 
  newLegend2 = newLegend
  rownames(newLegend2) = newLegend2$Factor = paste0("Spain_",sapply(rownames(newLegend),getPopnLabel,popnLabels)[1,])
  recieverLegend2 = rbind(recieverLegend2[!grepl("Spain",rownames(recieverLegend2)),],newLegend,newLegend2)
}

if(version_string=="Spain_A2_context_v7_v10a_contextPainting_v7_DONORpopn.v9"){
  
  popnLabels = popnLabels[!grepl("Spain",names(popnLabels))]
  popnLabelsNew = append(popnLabels,SPAINpopnLabels[grepl("cluster",names(SPAINpopnLabels))])
  names(popnLabelsNew) = gsub("cluster_","Spain_",names(popnLabelsNew))
    
  recieverLegendNEW = recieverLegendOrig[!grepl("Spain", recieverLegendOrig$Factor),]
  recieverLegendNEW = rbind(recieverLegendNEW,recieverLegend2[grepl('cluster',rownames(recieverLegend2)),])
  rownames(recieverLegendNEW) = recieverLegendNEW$Factor = gsub("cluster_","SpainPopn_",rownames(recieverLegendNEW))
  
  recieverLegend2 = recieverLegendNEW
  popnLabels = popnLabelsNew
}

if(version_string%in%c('Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12','Spain_A2_v1_v7a.new_contextPainting_v10_DONORpopn.v12.galCombine')){
    popnLabels = popnLabels[!grepl("Spain",names(popnLabels))]
    popnLabels$Basque = c(318,3,18)
    SPAINpopnLabels$Spain_Galicia_coast=unique(SPAINpopnLabels$Spain_Galicia_coast)
    popnLabelsNew = append(popnLabels,SPAINpopnLabels)
    
    # set up portugal labels
    popnLabelsNew$Spain_Portugal = 100
    recieverLegendNEW = rbind(DONORLegend,recieverLegend2)
    recieverLegend2 = rbind(recieverLegendNEW,makeFactorLegend("SpainPopn_100"))
    recieverLegend2["SpainPopn_100","colour"]="purple"
    recieverLegend2["SpainPopn_100","shape"]=25
    popnLabels = popnLabelsNew
}

######################################## 
# CHECK within-popn homogeneity
######################################## 

# compute within popn copy vector means
meanPopnCopyVectors <- matrix(NA,nrow=length(ALLpopns),ncol=dim(combinedChunkLengths)[2])
    rownames(meanPopnCopyVectors) <- ALLpopns
    colnames(meanPopnCopyVectors) <- colnames(combinedChunkLengths)
meanPopnCopydiffs = c()

for (Pop in ALLpopns){
  print(Pop)
  pop = as.numeric(str_split(Pop,"_",2)[[1]][2])
  if(length(grep("donor",Pop))==1) inds = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==pop]
  if(length(grep("Spain",Pop))==1) {
    inds = rownames(SPAIN.split.matrix)[SPAIN.split.matrix[,2]==pop]
    if(sum(grepl("SPAIN_A2",inds))==0) inds[!grepl("POPRES",inds)] = paste("SPAIN_A2_",inds[!grepl("POPRES",inds)],sep="")
  }

  meanPopnCopyVectors[Pop,] = colMeans(combinedChunkLengths[inds,])
  
  # get distance from the mean for all samples in popn==pop
  this = combinedChunkLengths[inds,]
  meanPopnCopydiffs = c(meanPopnCopydiffs,apply(this,MARGIN=1,FUN=norm_vec))
  
}


### save output ###

save(SUMMARYmat,SUMMARYmat_1,SUMMARYmat_1mean,SUMMARYmat_mean,SUMMARYmat_inds,mixmat,mixmat_output,mixmat_output_inds,mixmat_inds,
          DONOR.splitLevel,DONOR.split.matrix,SPAIN.split.matrix,SPAIN.splitLevel,summaryOutMeanChunkLengths,
          percentile05bootstrap,percentile95bootstrap,popnLabels,recieverLegend2,numToColour,SpainTree,meanPopnCopydiffs,meanPopnCopyVectors,
          combinedChunkLengths,combinedMeanChunkLengths,
     file=paste("../mixtureModel/mixtureModelOutput.",version_string,".Rdata",sep=""))

# if we have a new spain split matrix
if(version_stringStart==""){
    write.table(SPAIN.split.matrix,file=paste('mixtureModelOutput.SPAIN.split.matrix.',version_string,'.txt',sep=""),quote=F)
    SCP(paste('/mixtureModel/mixtureModelOutput.SPAIN.split.matrix.',version_string,'.txt',sep=""),direction="to",spain=T,recievedir="globetrotter_files/popnFiles")
    
    #write.table(SPAIN.split.matrix,file=paste('SPAIN.split.matrix_v',spain.chromv,'_v',spain.finev,'_',SPAIN.splitLevel2,extra,minSpainSplits,'.txt',sep=""),quote=F)
    SCP(paste('mixtureModel/SPAIN.split.matrix_v',spain.chromv,'_v',spain.finev,'_',SPAIN.splitLevel2,extra,minSpainSplits,'.txt',sep=""),direction="to",spain=T,recievedir="globetrotter_files/popnFiles")
    SCP(paste("mixtureModel/mixtureModelOutput.",version_string,".Rdata",sep=""),direction="to",spain=T,recievedir="globetrotter_files/popnFiles")
}
