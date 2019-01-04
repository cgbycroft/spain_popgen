#########################################
## BEGIN PROGRAM:

options(scipen=999)
library(nnls)

#######################
## SIMON'S FUNCTIONS:

getfit=function(predmat,fitdata,restrict=1){

  temp=predmat[-restrict,]
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]

  fitdata2=fitdata-predmat[restrict,]

  v=nnls(t(temp),fitdata2)
  x=v$x
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  v$x=newx
  #names(v$x)=colnames(predmat)
  names(v$x)=rownames(predmat)

  return(v)
}

getoverallfit=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  while(rep==1){
    q=getfit(predmat,fitdata,restrict=i)

    if(q$x[i]>0) rep=0
    i=i+1
  }

  return(q)
}

##############################################
## COMMAND TO RUN:

#ourmix=getoverallfit(donor.mat,recipient.vec)$x
#ourmix=ourmix[ourmix>0]
#ourmix=ourmix/sum(ourmix)

## END PROGRAM
###############################################
