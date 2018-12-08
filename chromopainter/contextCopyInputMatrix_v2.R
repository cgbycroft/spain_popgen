#########
# RScript to create input matrix for the context run of fineSTRUCTURE

library(dplyr)

args = commandArgs(trailingOnly=T)

# using inds as donors
#args = c(    "/well/donnelly/spain_structure/chromopainter_files/LOO-SP_PR_HM_NA/v7/combined2","../globetrotter_files/popnFiles/mixtureModelOutput.Spain_v1_v7a_contextPainting_v7_DONORpopn.v9.Rdata","SPAIN_A2_context/input/v2")

print(args)

load(args[2],verbose=T) # load list of donor and recipients to use

files = list.files(pattern="chunklengths.out",path=args[1],full.names=T)
files = files[grep("noselfcopy",files,invert=T)]

# combine the matrix with the others

lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

filesList = list()
rows = c()
for(file in files){
    print(file)
    a =   read.table(file,stringsAsFactors=F,header=T)
    rownames(a) = a[,1]
    a = a[,-1]        
    filesList = lappend(filesList,as.data.frame(a,stringsAsFactors=F))
    rows = c(rows,rownames(a))
}       
combinedMat = bind_rows(filesList)
rownames(combinedMat) = rows
combinedMat = as.matrix(combinedMat)

# note that for 28 donor samples the values will be 0 for all Spanish samples if we use the LOO versions. This is okay (I think).

# remove any NAs
notPort=!grepl("POPRES",rownames(SPAIN.split.matrix))
rownames(SPAIN.split.matrix)[notPort] = paste0("SPAIN_A2_",rownames(SPAIN.split.matrix)[notPort])

toRemoveDonor = unique(names(c(which(colSums(is.na(combinedMat))>0),
    which(!colnames(combinedMat)%in%rownames(DONOR.split.matrix)))))

combinedMat2 = combinedMat[!rownames(combinedMat)%in%toRemoveDonor,!colnames(combinedMat)%in%toRemoveDonor]

donors = rownames(combinedMat2)[rownames(combinedMat2)%in%rownames(DONOR.split.matrix)]
recip = rownames(combinedMat2)[rownames(combinedMat2)%in%rownames(SPAIN.split.matrix)]

# construct square matrix
combinedMat3 = matrix(0,nrow=length(donors)+length(recip),ncol=length(donors)+length(recip),)
rownames(combinedMat3) = colnames(combinedMat3) =c(recip, donors)

# populate Spanish parts (the rest is still zeros)
combinedMat3[recip,donors] = combinedMat2[recip,donors]



############## THIS IS WHERE THE MAGIC HAPPENS ##############
    
# v1: compute popn means based on donor popn groups

#fill = colMeans(combinedMat3[recip,donors])
#for(i in donors){
#    combinedMat3[i,donors[donors!=i]] = fill[donors!=i]
#}

# v2: Create vectors based on donor group assignment
dons = unique(DONOR.split.matrix[,2])
donorGroups = paste0("donorPopn_",dons)
ns = length(recip)
nt = dim(combinedMat3)[1]

for (d in dons){
    these = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==d]
    these = these[these%in%colnames(combinedMat3)] # exclude removed samples
    ni = length(these)
    spi = mean(combinedMat3[recip,these]) # take mean of spanish copying from all samples in donocr group
    zi = spi*(nt - ns - 1)/(ni - 1)

    combinedMat3[these,these] = zi
}

# set diagonal back to 0
diag(combinedMat3) = 0

# check for consistency in column means (between spain and donors)

a = mean(combinedMat3[recip,these])
b = mean(combinedMat3[donors,these])


#combinedMat3 = cbind(rownames(combinedMat3),combinedMat3)

# create force file information
# v1: create force file with continents: each with just one individual
#forcePopns = paste0("*",donors,"(",donors,")")

# v2: create force file with continents: one for each donor group
forcePopns = c()
for (d in dons){
    these = rownames(DONOR.split.matrix)[DONOR.split.matrix[,2]==d]
    these = these[these%in%colnames(combinedMat3)] # exclude removed samples
    
forcePopns = c(forcePopns,paste0("*donorPopn_",d,"(",paste(these,collapse=","),")"))

}

############## SAVE DATA  ##############

print("saving data")

print(paste0(args[3],"_forceFile.txt"))
print(paste0(args[3],"_coancestryMatrix.txt"))

# save force file
write.table(forcePopns,file = paste0(args[3],"_forceFile.txt"),quote=F,row.names=F,col.names=F)

# save the matrix
datafile <- file(paste0(args[3],"_coancestryMatrix.txt"), open = 'wt')    
on.exit(close(datafile))  # close on exit 
writeLines(c("Recipient",colnames(combinedMat3)),con=datafile, sep=' ')
writeLines('', con=datafile, sep='\n')
    
write.table(combinedMat3,datafile,quote=F,row.names=T,col.names=F)
        

