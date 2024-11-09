# Load libraries
# set libpath to find R libs in conda env
#.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/kosarek/miniconda3/envs/trax_env/lib/R/library",.libPaths()))
.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env/lib/R/library",.libPaths()))


library("DESeq")
args = commandArgs(trailingOnly = TRUE)
experimentname = args[1]
inputtable = args[2]
samplefile = args[3]

countTable = read.table(inputtable,check.names=FALSE)
conditions = colnames(countTable)

cds = newCountDataSet( countTable, conditions )
cds = estimateSizeFactors( cds )
scalingfactors = t(sizeFactors( cds ))
#colnames(scalingfactors) <- conditions
#write.table(rbind(colnames(readcounts),dds$sizeFactor),file=paste(experimentname,"/",experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)

write.table(rbind(conditions,scalingfactors),file=paste(experimentname,"/",experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)

normalizedrnas = sweep(countTable,2,scalingfactors, "/" )
write.table(normalizedrnas,paste(experimentname,"/",experimentname,"-normalizedreadcounts.txt", sep = ""), sep = "\t")
