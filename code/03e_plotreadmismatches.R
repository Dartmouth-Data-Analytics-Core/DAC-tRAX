# Load libraries
# set libpath to find R libs in conda env
#.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/kosarek/miniconda3/envs/trax_env/lib/R/library",.libPaths()))
.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env/lib/R/library",.libPaths()))

# Load libraries
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)

# Set command line arguments
args <- commandArgs(trailingOnly = TRUE)
outtype = ""
if(length(args) > 2){
outtype = args[3]
}

# Read in the data
counts <- read.table(args[1], header = TRUE, stringsAsFactors  = FALSE, row.names = NULL)

# Copy the data
selectcounts = counts


# Create a temp dataframe
temp = selectcounts

# Pivot the data from wide to long format
countsmelt = melt(temp, id.vars = c('count','type'))

# Aggregate counts
sampletotals = aggregate(countsmelt$value, list(countsmelt$variable), sum)

# Color palettes
colourCount = length(unique(countsmelt$seq))+1
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#unique(head(countsmelt))

#   other #f2dab2
#   tRNA #6666cc
#   pretRNA #a1ade5
#   miRNA #b4008d
#   snRNA #74531d
#   Mt_tRNA #87e275
#   Mt_rRNA #00754a
#   rRNA    #a5e5d9
#   snoRNA  #ff9e18
#   misc_RNA    #b2b2b2


#print(unique(countsmelt$seq))


# Print top of melted data
print(head(countsmelt))

# Plot visualization
ggplot(countsmelt,aes(x = variable, y = value,fill = count, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "fill",stat="identity") +
    geom_bar(position = "fill",stat="identity",color="black",show.legend=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    facet_wrap(~ type) +
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Percentage of Total Reads") + 
    labs(fill="Read\nMismatches")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    #scale_fill_manual(values = getPalette(colourCount))+
    #scale_fill_manual(values = typepal)+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

# Save visualization
ggsave(filename=args[2])
    
