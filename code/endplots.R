library(ggplot2)
library(reshape2)
library(scales)
library(plyr)
library(gridExtra)
library(getopt)

args <- commandArgs(trailingOnly = TRUE)



#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="tRNA.basePosition", value.name="read.density")


spec <- matrix(c(
        'cov'     , 'c', 1, "character", "coverage file from getcoverage.py (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'samples'    , 's', 1, "character", "sample file (required)",
        'uniquename'    , 'u', 1, "character", "file header for uniqifying files",
        'directory'    , 'f', 1, "character", "directory to place amino acid files",
        'allcov'    , 'a', 1, "character", "output coverages for all tRNAs (optional)",
        'multicov'    , 'm', 1, "character", "output coverages for all tRNAs on multiple pages(optional)",
        
        'combinecov'    , 'o', 1, "character", "output coverages for tRNAs combined",
        'mismatchcov'    , 'i', 1, "character", "output coverages for mismatches",
        'modomics'    , 'd', 2, "character", "optional modomics file",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);


coverages <- read.table(opt$cov, header = TRUE)
trnatable <- read.table(opt$trna)
sampletable <- read.table(opt$samples)
mismatchcov <- read.table(opt$mismatchcov)

modomicstable <- data.frame()



outputformat <- ".pdf"





#modomicstable <- read.table("/data/scratch/importrnaseq/agingtest/sacCer3-modomics.txt", header = TRUE)

outputfile <- opt$allcov
combinedfile <- opt$combinecov
multipage <- opt$multicov
#colnames(coverages)

myBreaks <- function(x){
    breaks <- c(min(0),floor(max(x)))
    names(breaks) <- attr(breaks,"labels")
    breaks
}
configurecov <- function(covplot){
covplot<- covplot+ theme_bw()
covplot<- covplot+theme(text = element_text(size = 4))
covplot<- covplot+theme(rect = element_rect(size=.01))
              
covplot<- covplot+theme(legend.key.size =  unit(.5, "cm")) # Change key size in the legend
covplot<- covplot+theme(legend.key.height =  unit(.5, "cm")) # Change key size in the legend
              
covplot<- covplot+theme(legend.text = element_text(size=6)) # Change the size labels in the legend
covplot<- covplot+theme(legend.title = element_text(size=6))             
              
covplot<- covplot+theme(axis.ticks = element_line(size = .1)) 
              
covplot<- covplot+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=3))

              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .1))
              
              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .15))
              
covplot<- covplot+theme(strip.switch.pad.grid = unit(.1, "cm"))
covplot<- covplot+theme(panel.margin = unit(.15, "lines"))
#This needs to be different from the combined version for some reason
covplot<- covplot+theme(axis.text.y=element_text(size=6))
            
covplot<- covplot+theme(strip.text.y = element_text(angle=0,size=10))
#as does this  
covplot<- covplot+theme(strip.text.x = element_text(angle=0,size=10))
covplot<- covplot+theme(axis.text=element_text(size=8))
covplot<- covplot+theme(axis.title=element_text(size=8))
#covplot<- covplot+scale_y_continuous(limits = c(0, 1))


covplot
}




#remove columns with too many NAs
coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8]
coveragemelt = melt(coverages, id.vars = c("Feature", "Sample"))
coveragemelt[is.na(coveragemelt)] <- 0
coveragemeltagg <- aggregate(coveragemelt$value, by=list(Feature = coveragemelt$Feature, Sample = sampletable[match(coveragemelt$Sample,sampletable[,1]),2], variable = coveragemelt$variable), FUN=mean)
colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
coveragemelt <- coveragemeltagg

write.table(coveragemelt,"aggtables.txt" ,sep = "\t")

acceptorType = trnatable[match(coveragemelt$Feature, trnatable[,1]),3]
acceptorType <- factor(acceptorType, levels = sort(unique(acceptorType)))


sortcovmelt <- coveragemelt[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature),]  
sortacceptor <- acceptorType[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature)]





modomicstable <- data.frame(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)

modomicstable <- read.table(text ="",col.names = c("trna", "mod", "pos"),colClasses = c("character", "character", "character")) #(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)
#colNames(modomicstable) <- c("trna", "mod","pos")
if(!is.null(opt$modomics) && file.exists(opt$modomics) ){
    modomicstable <- read.table(opt$modomics, header = TRUE)
    modomicstable$pos <- paste("X",modomicstable$pos, sep = "")
}



modomicstable <- modomicstable[as.character(modomicstable$pos) %in% unique(as.character(coveragemelt$variable)),]
modomicstable$dist <- match(modomicstable$pos,  levels(coveragemelt$variable)) #factor(modomicstable$pos, levels = levels(coveragemelt$variable))
modomicstable$Feature <- factor(modomicstable$trna,levels = levels(coveragemelt$Feature) )

stopmods = c("m1A","m2,2G","m1G","m1I","m3C")
modomicstable <- modomicstable[modomicstable$mod %in% stopmods,]

colnames(modomicstable)[colnames(modomicstable) == 'mod'] <- 'Modification'


    
    
#modomicstable$dist <- modomicstable$pos
#head(coveragemelt)
#unique(as.character(coveragemelt$variable))
#
#geom_bar(aes(fill = factor(baseMod, levels=c("m1A", "m1G", "m3C", "other bases", "not documented"))), stat="identity")
#modomicstable <- data.frame(Feature = factor(c("nmt-tRNA-Leu-TAA-1"),levels(coveragemelt$Feature)), pos =10)
#ggsave(filename=outputfile, width = 30, height = 30)
#3*length(unique(coveragemelt$Feature))
#5*length(unique(coveragemelt$Sample))
#dev.off()
scalefactor = .5

#ggsave(filename=outputfile, allcoverages,height=scalefactor*1*length(unique(coveragemelt$Feature)),width=2+scalefactor*5*length(unique(coveragemelt$Sample)), limitsize=FALSE, dpi = 600)
#ggsave(filename=outputfile, allcoverages,dpi = 1000, limitsize=FALSE)
#set dpi

#unique(acceptorType)
if(!is.null(opt$directory) && FALSE){
for (curramino in unique(acceptorType)){
#

aminodata = coveragemelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]
aminocoverage <- ggplot(aminodata,aes(x=variable,y=value), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE)+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+ ylab("Percentage of Read Ends") + xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c("Start","m1g","m22g","anticodon","varloop","m1a","65","tail")) + guide_legend(title="RNA modifications")   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))

aminoname = paste(opt$directory,"/",curramino,"-endcoverage",outputformat, sep = "")
#print(aminoname)
ggsave(filename=aminoname, aminocoverage,height=2 + scalefactor*1.5*length(unique(aminodata$Feature)),width=scalefactor*5*length(unique(aminodata$Sample)), limitsize=FALSE, dpi = 600)

}
}

if(!is.null(opt$mismatchcov))
{
mismatchcovs <- read.table(paste(opt$mismatchcov), header = TRUE)




mismatchcovs <- mismatchcovs[ , colSums(is.na(mismatchcovs)) < nrow(mismatchcovs)/8]

mismatchmelt = melt(mismatchcovs, id.vars = c("Feature", "Sample"))
mismatchmelt[is.na(mismatchmelt)] <- 0


mismatchmeltagg <- aggregate(mismatchmelt$value, by=list(Feature = mismatchmelt$Feature, Sample = sampletable[match(mismatchmelt$Sample,sampletable[,1]),2],variable = mismatchmelt$variable), FUN=mean)
colnames(mismatchmeltagg)[colnames(mismatchmeltagg) == "x"]  <- "value"
mismatchmeltagg$Sample <- factor(mismatchmeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
mismatchmelt <- mismatchmeltagg



mismatchcoverage <- ggplot(mismatchmelt,aes(x=variable,y=value, size = 2)) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = 8,angle=0)) + ylab("Percent Mismatches") +   xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) + labs(fill="Mappability", vline="RNA\nModification") 




smallcovall <- mismatchcoverage
smallcovall <- configurecov(mismatchcoverage)
ggsave(filename=paste(opt$unique, "-mismatches.pdf",sep= ""), smallcovall,height=scalefactor*(2 + 1.5*length(unique(mismatchmelt$Feature))),width=scalefactor*(2+5*length(unique(mismatchmelt$Sample))), limitsize=FALSE, dpi = 600)




for (curramino in unique(acceptorType)){

aminodata = mismatchmelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]


aminocoverage   <- ggplot(aminodata,aes(x=variable,y=value), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Percent Mismatches") + xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X13","X22","X31","X39","X53","X61","X73"), labels=c("Start","D-loop start","D-loop end","AC-loop start","AC-loop end","T-loop start","T-loop end","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))

aminocoverage <- configurecov(aminocoverage)

aminoname = paste(opt$directory,"/",curramino,"-mismatches",outputformat, sep = "")
#print(aminoname)
ggsave(filename=aminoname, aminocoverage,height=scalefactor*(2 + 1.5*length(unique(aminodata$Feature))),width=scalefactor*(2+5*length(unique(aminodata$Sample))), limitsize=FALSE, dpi = 600)


}
}

if(!is.null(opt$unique)){


multamino <- read.table(paste(opt$unique, "-multaminoends.txt",sep= ""), header = TRUE)
multactable <- read.table(paste(opt$unique, "-multacends.txt",sep= ""), header = TRUE)
multtrnas <- read.table(paste(opt$unique, "-multtrnaends.txt",sep= ""), header = TRUE)
uniquetable <- read.table(paste(opt$unique, "-uniqueends.txt",sep= ""), header = TRUE)

uniquetable$maptype = "Transcript specific"
multactable$maptype = "Isotype Specific"
multamino$maptype = "Nonspecific"
multtrnas$maptype = "Isodecoder Specific" 
#allmulttables <- rbind(uniquetable,multactable,multamino,multtrnas)

#allmulttables <- rbind(multamino,multtrnas,multactable,uniquetable)
allmulttables <- rbind(uniquetable, multtrnas,multactable,multamino)


#print(head(allmulttables))
allmulttables <- allmulttables[ , colSums(is.na(allmulttables)) < nrow(allmulttables)/8]
#allmulttables <- allmulttables[ , colSums(is.na(allmulttables)) < nrow(allmulttables)/8]
#print(colnames(allmulttables))
allmultmelt = melt(allmulttables, id.vars = c("Feature", "Sample", "maptype"))
allmultmelt[is.na(allmultmelt)] <- 0


allmultmeltagg <- aggregate(allmultmelt$value, by=list(Feature = allmultmelt$Feature, Sample = sampletable[match(allmultmelt$Sample,sampletable[,1]),2], maptype = allmultmelt$maptype,variable = allmultmelt$variable), FUN=mean)
colnames(allmultmeltagg)[colnames(allmultmeltagg) == "x"]  <- "value"
allmultmeltagg$Sample <- factor(allmultmeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
allmultmelt <- allmultmeltagg

#
#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Unique tRNA","Unique Decoder","Unique Acceptor","Nonunique Acceptor"))
#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Nonunique Acceptor","Unique Acceptor","Unique Decoder","Unique tRNA"))
#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Transcript specific","Isodecoder Specific","Isotype Specific","Nonspecific"))
allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Nonspecific","Isotype Specific","Isodecoder Specific","Transcript specific"))


allmultmelt <- allmultmelt[order(-as.numeric(allmultmelt$maptype)),]
#print(head(allmultmelt))
  



for (curramino in unique(acceptorType)){

aminodata = allmultmelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]
#This bit messes up the sample ordering
#aminodata$Sample <- gsub("_", " ", aminodata$Sample)  



aminocoverage   <- ggplot(aminodata,aes(x=variable,y=value, fill = maptype, order=-as.numeric(maptype)), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Percent Mismatches") + xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X13","X22","X31","X39","X53","X61","X73"), labels=c("Start","D-loop start","D-loop end","AC-loop start","AC-loop end","T-loop start","T-loop end","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))
#all formula


aminocoverage <- configurecov(aminocoverage)

aminoname = paste(opt$directory,"/",curramino,"-uniqends",outputformat, sep = "")
#print(aminoname)
ggsave(filename=aminoname, aminocoverage,height=scalefactor*(2 + 1.5*length(unique(aminodata$Feature))),width=scalefactor*(2+5*length(unique(aminodata$Sample))), limitsize=FALSE, dpi = 600)


}



uniqcoverage <- ggplot(allmultmelt,aes(x=variable,y=value, fill = maptype, order=-as.numeric(maptype)), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = 8,angle=0)) + ylab("Percent Mismatches") +   xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ guide_legend(title="RNA modifications")   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))




smallcovall<-uniqcoverage
configurecov <- configurecov(smallcovall)
ggsave(filename=paste(opt$unique, "-mismatches2.pdf",sep= ""), uniqcoverage,height=scalefactor*(2 + 1.5*length(unique(allmultmelt$Feature))),width=scalefactor*(2+5*length(unique(allmultmelt$Sample))), limitsize=FALSE, dpi = 600)


}




if(!is.null(opt$indelcov))
{
indelcovs <- read.table(paste(opt$indelcov), header = TRUE)




indelcovs <- indelcovs[ , colSums(is.na(mismatchcovs)) < nrow(mismatchcovs)/8]

indelmelt = melt(indelcovs, id.vars = c("Feature", "Sample"))
indelmelt[is.na(indelmelt)] <- 0


indelmeltagg <- aggregate(indelmelt$value, by=list(Feature = indelmelt$Feature, Sample = sampletable[match(indelmelt$Sample,sampletable[,1]),2],variable = indelmelt$variable), FUN=mean)
colnames(indelmeltagg)[colnames(indelmeltagg) == "x"]  <- "value"
indelmeltagg$Sample <- factor(indelmeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
indelmelt <- indelmeltagg



indelcoverage <- ggplot(indelmelt,aes(x=variable,y=value, size = 2)) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = modomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = 8,angle=0)) + ylab("Percent Deletions") +   xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) + labs(fill="Mappability", vline="RNA\nModification") 




smallindelall <- configurecov(indelcoverage)
ggsave(filename=paste(opt$unique, "-indels.pdf",sep= ""), smallindelall,height=scalefactor*(2 + 1.5*length(unique(mismatchmelt$Feature))),width=scalefactor*(2+5*length(unique(mismatchmelt$Sample))), limitsize=FALSE, dpi = 600)




for (curramino in unique(acceptorType)){

aminodata = indelmelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]


aminocoverage   <- ggplot(aminodata,aes(x=variable,y=value), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Percent Mismatches") + xlab("tRNA position") + scale_y_continuous(limits = c(0, 1),breaks=c(0,1),labels=c("0%","100%")) +scale_x_discrete(breaks=c("X1","X13","X22","X31","X39","X53","X61","X73"), labels=c("Start","D-loop start","D-loop end","AC-loop start","AC-loop end","T-loop start","T-loop end","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon", "tail"))

aminocoverage <- configurecov(aminocoverage)

aminoname = paste(opt$directory,"/",curramino,"-indels",outputformat, sep = "")
#print(aminoname)
ggsave(filename=aminoname, aminocoverage,height=scalefactor*(2 + 1.5*length(unique(aminodata$Feature))),width=scalefactor*(2+5*length(unique(aminodata$Sample))), limitsize=FALSE, dpi = 600)


}



}





