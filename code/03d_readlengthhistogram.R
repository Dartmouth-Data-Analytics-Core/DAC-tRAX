# Load libraries
# set libpath to find R libs in conda env
#.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/kosarek/miniconda3/envs/trax_env/lib/R/library",.libPaths()))
.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env/lib/R/library",.libPaths()))

# Load librarues
library(ggplot2)
library(reshape2)
library(scales)

# Set command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Read in read lengths and sample table
readlengths <- read.table(args[1], header = TRUE)
sampledata <- read.table(args[2])

# Plot visualization
if(FALSE){
ggplot(readlengths,aes(x = variable, y = value,fill = seq, stat="identity")) +
	geom_bar(position = "fill",stat="identity") + 
    geom_bar(position = "fill",stat="identity",color="black",show_guide=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    xlab("Sample") +
    ylab("Percentage of Reads") + 
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5))
    }
    

	

binw <- 2
maxlen <- min(c(200, max(readlengths[,1])))
readlengths$bin <- floor(readlengths$Length / binw)*binw

# Aggregate binned read lengths for tRNAs
binnedtrnareadlengths <- aggregate(trnas ~ bin * Sample, readlengths, sum)

# Aggregate binned read lengths for non tRNAs
binnedreadlengths <- aggregate(other ~ bin * Sample ,readlengths, sum)

# Merge
temp = merge(binnedtrnareadlengths,binnedreadlengths,by=c("bin", "Sample"))

# Pivot from wide to long format
binnedlengths = melt(temp, id.vars = c('bin', 'Sample'))

# Filter and order
binnedlengths[is.na(binnedlengths)] <- 0
binnedlengths <- binnedlengths[order(binnedlengths$Sample,binnedlengths$bin),]

# Factor
binnedlengths$variable <- factor(binnedlengths$variable, levels = c("other", "trnas"))
binnedlengths$Sample <- factor(binnedlengths$Sample, levels = sampledata[,1])

# Plot visualization
ggplot(data=binnedlengths, aes(x = bin, y = value, fill = variable))+
    geom_bar(stat="identity")+
    facet_wrap( ~ Sample, scales="free",ncol = 3)+ 
    xlim(0,maxlen) + 
    ylab("Count") + 
    xlab("Read Length") + 
    scale_fill_discrete(name="Gene Type") + 
    scale_x_continuous(labels = comma, limits=c(0, maxlen)) + 
    theme_bw()

# Save visualization
ggsave(filename=args[3],limitsize=FALSE,width = 8, height = .5 * length(unique(binnedlengths$Sample)))#, width = 3 * length(unique(binnedlengths$Sample))) 
    
