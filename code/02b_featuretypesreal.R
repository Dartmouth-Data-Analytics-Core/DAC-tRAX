# Load libraries
# set libpath to find R libs in conda env
#.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/kosarek/miniconda3/envs/trax_env/lib/R/library",.libPaths()))
#.libPaths(c("/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env/lib/R/library",.libPaths()))

#----- Loading libraries
library(ggplot2)
library(reshape2)
library(scales)

# Retrieve the command-line arguments passed to the Rscript
args <- commandArgs(trailingOnly = TRUE)

# Set output type initialization
outtype = ""
if(length(args) > 2){
outtype = args[3]
}


#Rscript trailerbarplot.R hg19-trailertable.txt hg19-barplot.png

#----- Read a table from the first command-line argument
counts <- read.table(args[1],check.names=FALSE)

# Create a copy of the counts dataframe to work with
selectcounts = counts

#Add a new column 'seq' to 'selectcounts' which contains row names as factors, ordered in reverse
temp = cbind(selectcounts, seq = factor(rownames(selectcounts),rev(rownames(selectcounts))))

# Uncommenting this line would reverse the factor order
#levels(temp$seq) <- rev(rownames(selectcounts))

# Melt the dataframe 'temp' ro convert it from wide format to long format
# keep 'seq' as an identifier column
countsmelt = melt(temp, id.vars = c('seq'))

# Update the factor levels of 'seq' in the 'countsmelt' dataframe to match the rownames of 'selectcounts'
countsmelt = within(countsmelt, seq <- factor(seq, 
    rownames(selectcounts)))

#head(countsmelt)

# Aggregate the values in 'countsmelt' by summing the 'value' column grouped by the 'variable' colummn
# The result is stored in 'sampletotals', which contains he total counts for each variable (sample)
sampletotals = aggregate(countsmelt$value, list(countsmelt$variable), sum)

# Filter the 'countsmelt' df to keep only rows where 'value' is greater than 0
countsmelt = countsmelt[countsmelt$value > 0,]

# Check if the output type is not set to 'all'
if(outtype != "all"){

# Further filter 'countsmelt' to keep only rows where 'value' is greater than 1% of the total count for the corresponding sample
countsmelt = countsmelt[countsmelt$value > sampletotals$x[countsmelt$variable] / 100,] #filters out those types below a certain level
}

# Plot visualizations
ggplot(countsmelt,aes(x = variable, y = value,fill = seq, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(stat="identity") + 
    geom_bar(stat="identity",color="black",show.legend=FALSE) + 
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Total Reads") + 
    labs(fill="Read\nType") + 
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5))

# Save visualization
ggsave(filename=args[2])
    
