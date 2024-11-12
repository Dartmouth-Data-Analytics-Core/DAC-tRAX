library(ggplot2)
library(reshape2)
library(scales)
library(getopt)


spec <- matrix(c(
        'cov'     , 'c', 1, "character", "coverage file from getcoverage.py (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'samples'    , 's', 1, "character", "sample file (required)",
        'allcov'    , 'a', 1, "character", "output coverages for all tRNAs (optional)",
        'multicov'    , 'm', 1, "character", "output coverages for all tRNAs on multiple pages(optional)",
        'directory'    , 'f', 1, "character", "directory to place amino acid files",        
        'combinecov'    , 'o', 1, "character", "output coverages for tRNAs combined",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (length(readLines(opt$cov, n = 2)) < 2){quit('yes')}

coverages <- read.table(opt$cov, header = TRUE)

trnatable <- read.table(opt$trna)
sampletable <- read.table(opt$samples)

outputfile <- opt$allcov
combinedfile <- opt$combinecov
multipage <- opt$multicov
#colnames(coverages)

expand.delimited <- function(x, col1=1, col2=2, sep=",") {
  rnum <- 1
  expand_row <- function(y) {
    factr <- y[col1]
    strng <- toString(y[col2])
    expand <- strsplit(strng, sep)[[1]]
    num <- length(expand)
    factor <- rep(factr,num)
    return(as.data.frame(cbind(factor,expand),
          row.names=seq(rnum:(rnum+num)-1)))
    rnum <- (rnum+num)-1
  }
  expanded <- apply(x,1,expand_row)
  df <- do.call("rbind", expanded)
  names(df) <- c(names(x)[col1],names(x)[col2])
  return(df)
}


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
#smallcovall< covplot+theme(axis.ticks.y=element_blank())
              
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

covplot
}
outputformat <- ".pdf"
outputindivformat <- ".png"
scalefactor = .5
coverages <- coverages[ , colSums(is.na(coverages)) < nrow(coverages)/8]


coveragemelt = melt(coverages, id.vars = c("Feature", "Sample"))


coveragemelt[is.na(coveragemelt)] <- 0


#Normalization
#This normalization takes way too long



coveragemeltagg <- aggregate(coveragemelt$value, by=list(Feature = coveragemelt$Feature, Sample = sampletable[match(coveragemelt$Sample,sampletable[,1]),2], variable = coveragemelt$variable), FUN=mean)
colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
coveragemelt <- coveragemeltagg




#vector of amino acids
locustable = expand.delimited(trnatable,3,2)
acceptorType = locustable[match(coveragemelt$Feature, locustable[,2]),1]

sortcovmelt <- coveragemelt[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature),]  
sortacceptor <- acceptorType[order(coveragemelt$variable, coveragemelt$Sample,as.numeric(acceptorType),coveragemelt$Feature)]



covsummary <- ggplot(sortcovmelt,aes(x=variable,y=value, fill = sortacceptor, order = as.numeric(sortacceptor))) + facet_grid( ~ Sample, scales="free") +xlab("Position")+ geom_bar(stat="identity")+theme(axis.text.y=element_text(colour="black",size=8), strip.text.y = element_text(angle=0,size=4),strip.text.x = element_text(angle=0,size=8),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8))+ ylab("Read Share") +   scale_y_continuous(breaks=myBreaks, labels = c("0","1")) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","tail")) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))
ggsave(filename=combinedfile,covsummary)

allmultmelt <- sortcovmelt
for (currtranscript in unique(allmultmelt$Feature)){

transcriptdata = allmultmelt[allmultmelt$Feature == currtranscript,]
 


transcriptcoverage   <- ggplot(transcriptdata,aes(x=variable,y=value), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X13","X22","X31","X39","X53","X61","X73"), labels=c("Start","D-loop start","D-loop end","acc-loop start","acc-loop end","D-loop start","D-loop end","tail")) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")
#all formula


transcriptcoverage <- configurecov(transcriptcoverage)
transcriptcoverage <- transcriptcoverage + theme(legend.box = "horizontal")

transcriptname = paste(opt$directory,"/indiv/",currtranscript,"-coverage",outputindivformat, sep = "")
#print(aminoname)
ggsave(filename=transcriptname, transcriptcoverage,height=scalefactor*(2 + 1.5*length(unique(transcriptdata$Feature))),width=scalefactor*(2+5*length(unique(transcriptdata$Sample))), limitsize=FALSE, dpi = 600)


}


locusplot <- ggplot(coveragemelt,aes(x=variable,y=value)) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=4),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0,size=6),strip.text.x = element_text(angle=0,size=8))+ ylab("read count") +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c("Start","anticodon","End"))  
ggsave(locusplot, filename=outputfile,height=.5*length(unique(coveragemelt$Feature)),width=2*length(unique(coveragemelt$Sample)), limitsize=FALSE)


