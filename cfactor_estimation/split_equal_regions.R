
# program to split chromosomes into contiguous regions of equal length in cM
setwd('/well/donnelly/spain_structure/chromopainter_files')
args <- commandArgs(trailingOnly = TRUE)

len = as.numeric(args[1])
snpsFiles= args[2]
outputDIR = args[3]
#snpsFiles= "SPAIN_A2/v1/chromofiles/SPAIN.A2.chr%%_v1.chromo.haps"
print(paste("Finding positions such that chromosomes are split into segments of length ",len,"cM.",sep=""))

#print(args)

for (chrom in c(1:22)){
    print(chrom)
    recombmap=read.table(paste('../phasing/ALL_1000G_phase1integrated_v3_impute_macGT1/genetic_map_chr',chrom,'_combined_b37_forchrompainter.txt',sep=""),header=T)
    print(head(recombmap))
    snpsFile=gsub("%%",chrom,snpsFiles)
    
    snps=read.table(snpsFile,header=F)
    #print(len)
    #print(snps[1:10,])
    maxlength = max(recombmap$Genetic_Map.cM.)
    #print(maxlength)
    #print(round(maxlength))
    segments = seq(0,round(maxlength),len)
    #print(segments)

    i=1
    index <- c()
    for (s in segments){
        print(s)
        snp <- which(abs(recombmap$Genetic_Map.cM.-s)==min(abs(recombmap$Genetic_Map.cM.-s)))
        index[i] <- snp[1]
        
        i=i+1
    }

    positions <- recombmap[index,"position"]

    position.pairs <- as.data.frame(t(c(positions[1],positions[2])))
    if (length(positions)>2){
        for (i in c(2:(length(positions)-1))){
            position.pairs <- rbind(position.pairs,c(positions[i]+1, positions[i+1]))
        }
    }
    
    # get snps numbers that match these positions
    position.pairs$snps1 <- NA
    position.pairs$snps2 <- NA
    for (p in c(1:(length(positions)-1))){
        pos <- position.pairs[p,1:2]
        these <- which((snps>=pos[,1])&(snps<=pos[,2]))
        position.pairs[p,3:4] <- c(min(these),max(these))
    }
    
    write.table(position.pairs,file=paste(outputDIR,'/chr',chrom,'_split-',len,'cM_position_pairs.txt',sep=""),quote=F,row.names=F,col.names=F)

}


print("Done!")
print(paste("Output files are in",outputDIR))
