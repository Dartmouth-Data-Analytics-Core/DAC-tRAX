# Load libraries
# set libpath to find R libs in conda env
#.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/kosarek/miniconda3/envs/trax_env/lib/R/library",.libPaths()))
.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env/lib/R/library",.libPaths()))

# Load Libraries
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)

# Set command lind arguments
args <- commandArgs(trailingOnly = TRUE)

# Read in the counts
counts <- read.table(args[1],check.names=FALSE)

# Copy the counts to a new variable
selectcounts = counts

# Create a temp
temp = cbind(selectcounts, seq = factor(rownames(selectcounts), rev(rownames(selectcounts)), ordered = TRUE))

# Collapse from wide format into long
countsmelt = melt(temp, id.vars = c('seq'))

# Factor the data
countsmelt = within(countsmelt, seq <- factor(seq, 
    rownames(selectcounts), 
    ordered = FALSE))

# Aggregate
sampletotals = aggregate(countsmelt$value, list(countsmelt$variable), sum)

# Filter the data
countsmelt = countsmelt[countsmelt$value > sampletotals$x[countsmelt$variable] / 100,]

# Set color palettes
colourCount = length(unique(countsmelt$seq))+1
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Set type color palette for plotting
typepal <- c(
  "tRNA" = "#6666cc", 
  "pretRNA" = "#a1ade5", 
  "miRNA" = "#b4008d",
  "snRNA" = "#74531d",
  "Mt_tRNA" = "#87e275",
  "Mt_rRNA" = "#00754a",
  "rRNA" = "#a5e5d9",
  "snoRNA" = "#ff9e18",
  "misc_RNA" = "#b2b2b2",
  "other" = "#f2dab2"
)

# Extract extra genes
extragenes = setdiff(unique(countsmelt$seq),names(typepal))

# Set color palette for "other"
otherpal = getPalette(length(extragenes))
names(otherpal) = extragenes

# Finalize the color palette    
typepal = c(typepal, otherpal)

# Plot visualization
ggplot(countsmelt,aes(x = variable, y = value,fill = seq, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "fill",stat="identity") +
    geom_bar(position = "fill",stat="identity",color="black",show.legend=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Percentage of Total Reads") + 
    labs(fill="Read\nType")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    #scale_fill_manual(values = getPalette(colourCount))+
    scale_fill_manual(values = typepal)+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

# Save visualization
ggsave(filename=args[2])
    
