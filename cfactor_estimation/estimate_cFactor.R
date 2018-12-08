# This script is used to estimate the c-factor for fineSTRUCTURE input when using chunklengths as the coancestry matrix.
# Author: Clare Bycroft
# Date: August 2014
##############
# Input:
# A series of chromopainter (combined) summary matrices of the form:
# <file_prefix>_chr14_segment_2.chunklengths.out
##############
# Output:
# - A single scalar value of the c-factor
# - A matrix of inflation factors for individuals i and j (c-factor is the average of these)
# NOTE: Output directory is the same as the directory for the input files.
##############
# Bash script arguments:
# Rscript <file_prefix (including input dir)>
##############

library(stringr)

############## Useful Functions ##############

fixNames <- function(x) {
        rownames(x) <- sub("X","",rownames(x))
        colnames(x) <- sub("X","",colnames(x))
        return(x)
      }

############## INPUT DATA ##############
args <- commandArgs(trailingOnly = TRUE)
#args="SPAIN_A2/v1_splits/combined/renamed/SPAIN.A2_v1_splits"
#args = "LOO-SP_PR_HM_NA/v7_splits/combined2/LOO-SP_PR_HM_NA_v7_splits"

print(args)
print("Reading in data...")

inputFilePrefix <- args[1]

#do we double the number of segments to account for diploidy?
if (args[2] == "TRUE") double = T
if (args[2] == "FALSE") double = F

print(paste("double segments? = ",double))

system(paste("ls ",inputFilePrefix,"*.chunklengths.out > segmentFiles.tmp",sep=""))

seg.Files <- read.table('segmentFiles.tmp',header=F,stringsAsFactors=F)

seg = 1

# how many chromosomes available? Might not be all 22!

strings= unique(str_extract(seg.Files[,1],"chr.*_"))
Chroms = str_extract(strings[!is.na(strings)],"\\d+")
    
for (chr in Chroms){

    files = seg.Files[grep(paste("chr",chr,"_",sep=""),seg.Files[,1]),]
    
    for (file in files){
        print(file)
        f = as.matrix(read.table(file,stringsAsFactors=F,header=T,row.names=1))
        #c <- fixNames(f)
        if (seg==1){
            nsamples = nrow(f)
            Data <- array(NA,c(nrow(seg.Files),nsamples,ncol(f)))
            NamesRows = rownames(fixNames(f))
            NamesCols = colnames(fixNames(f))
        }
       # print(head(f))
       # print(seg)
        Data[seg,,] <- f
        seg = seg + 1
    }
}




############## Estimating c-factor ##############

print("Calculating theoretical and empirical variance...")

Sij = t(sapply(c(1:dim(Data)[2]),FUN=function(i) {
    
    d = Data[,i,] # get vectors across segments for sample x
    Si = colSums(d) # sum ind i with ind j over all segments
}))

Sij2 = t(sapply(c(1:dim(Data)[2]),FUN=function(i) {

    d = Data[,i,] # get vectors across segments for sample x
    Si2 = colSums(d*d) # sum over all segments, the square of ind i with ind j, segment k
}))

Pij = Sij/rowSums(Sij)  # note that rowSums(Sij) is the same (approx) when using chunklengths, as everybody copies from all the genome.

Ri = dim(Data)[1] # number of segments per individual i. In this case with chunklengths it is all the same (x 2 because two sets of chromosomes).  rowSums(Sij)/Ri is just the length of each segment (defined when breaking up chromosomes) but on average Ri is slightly shorter. e.g if segments set to 40cM, rowSums(Sij)/Ri = 39.6.

if (double==T) Ri = dim(Data)[1]*2

#T = rowSums(Sij)[1]
    
N = dim(Data)[2]

empiricalVar = Sij2/(Ri-1) - (Sij*Sij)/(Ri*(Ri-1))
theoreticalVar = rowSums(Sij)*Pij*(1-Pij)/Ri
theoreticalVar2 = Sij*(1-Pij)/Ri

print("Calculating ratio of the two to get c-factor...")

Cij = 2*empiricalVar/theoreticalVar # Factor of 2 only necessary when matrix is square. i.e not in panel mode. ==> halve the number if this is the case.


c = sum(rowSums(Cij,na.rm=T))*(1/(N*(N-1)))
print(paste("c-factor estimate is ",c))

############### Write data to files ##############

rownames(Cij) = NamesRows
colnames(Cij) = NamesCols

if (double==F) {
    cFactorFile = paste(inputFilePrefix,"_chunklengths.c-factor.txt",sep="")
    CijFile =  paste(inputFilePrefix,"_chunklengths.variance.inflation.out",sep="")
}
if (double==T) {
    cFactorFile = paste(inputFilePrefix,"_chunklengths.c-factor2.txt",sep="")
    CijFile =  paste(inputFilePrefix,"_chunklengths.variance.inflation2.out",sep="")
}

write.table(c,file = cFactorFile,quote=F,col.names=F,row.names=F)
write.table(Cij,file = CijFile,quote=F)

print(paste("c-factor value is saved in ",cFactorFile,sep=""))
print(paste("Pairwise variance inflation matrix is saved in ",CijFile,sep=""))

system("rm segmentFiles.tmp")

print("DONE!")
print("Stats is cool.")      
