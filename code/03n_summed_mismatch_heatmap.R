# Load necessary libraries
library(reshape2)
library(getopt)
library(ggplot2)

# Set up command line arguments
args <- commandArgs(trailingOnly = TRUE)
spec <- matrix(c(
  'runname', 'n', 1, "character", "name of exp",
  'mismatch', 'm', 1, "character", "mismatches",
  "pairfile", 'p', 1, "character", "pairfile",
  'directory', 'd', 1, "character", "output directory (required)"
), ncol = 5, byrow = TRUE)

opt = getopt(spec)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in the coverage data and factor positions
data <- read.table(opt$mismatch, header = TRUE)
data$position <- factor(data$position, levels = unique(data$position))

# Read in the replicate file and set column names
repfile <- read.csv(opt$pairfile, sep = " ", header = FALSE)
colnames(repfile) <- c("SampleID", "Group", "Fastq")

# Process data by unique groups
uniqueGroups <- unique(repfile$Group)
uniqueSamples <- list()
for(i in uniqueGroups) {
  subset <- unique(repfile[repfile$Group == i,]$SampleID)
  uniqueSamples[[i]] <- subset
}

# Split data by sample
df_subsets <- split(data, data$Sample)
keep <- c("Feature", "Sample", "position", "mismatchedbases")

# Transform data to wide format
processDFs <- function(df) {
  df <- df[, colnames(df) %in% keep]
  wide_df <- dcast(df, Feature ~ position, value.var = "mismatchedbases", fill = 0)
  rownames(wide_df) <- wide_df$Feature
  wide_df$Feature <- NULL
  return(wide_df)
}
wide_df_list <- lapply(df_subsets, processDFs)

# Summing groups
summedGroups <- list()
for (i in names(uniqueSamples)) {
  sampleIDs <- uniqueSamples[[i]]
  groupDFs <- wide_df_list[sampleIDs]
  for (j in names(groupDFs)) {
    groupDFs[[j]] <- groupDFs[[j]][, !colnames(groupDFs[[j]]) %in% "Sample"]
  }
  summedDF <- Reduce('+', groupDFs)
  summedGroups[[i]] <- summedDF
}

# Calculate the difference between groups
groupNames <- names(uniqueSamples)
group1 <- groupNames[[1]]
group2 <- groupNames[[2]]
result <- as.data.frame(summedGroups[[group1]] - summedGroups[[group2]])

# Clean up results
rownames(result) <- gsub("tRNA-", "", rownames(result))
aminoAcids <- as.data.frame(sub("-.*", "", rownames(result)))
colnames(aminoAcids) <- "Amino_Acid_Family"

# Prepare data for heatmap
heatmapData <- as.matrix(result)
scaledData <- as.data.frame(t(scale(t(heatmapData))))
scaledData$Iso <- rownames(scaledData)
scaledData <- melt(scaledData)
scaledData$Acceptor <- sub("-.*", "", scaledData$Iso)

# Color palette
col <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
         "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
         "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC",
         "#E5D8BD", "#FDDAEC", "#F2F2F2", "#1B9E77", "#D95F02", "#7570B3", "#E7298A")
unique_aa <- unique(aminoAcids$Amino_Acid_Family)
colors <- setNames(col[seq_along(unique_aa)], unique_aa)

# Plot heatmap
directory <- opt$directory
expname <- opt$runname
posname <- paste(directory, "/", expname, "-Full_relative_summed_tRNA_mismatches.png", sep = "")

hm <- ggplot(scaledData, aes(x = variable, y = Iso, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#6A9FB5", mid = "#FFFFF9", high = "#D46A6A", midpoint = 0) +
  facet_grid(Acceptor ~ ., scales = "free") +
  theme_classic() +
  labs(y = "", x = "Position", fill = "Z-Score") +
  theme(axis.text.y = element_text(size = 6, face = "bold"),
        axis.title.x = element_text(size = 26, face = "bold"),
        strip.text = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size = 16, face = "bold", angle = 90),
        legend.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        title = element_text(size = 14, face = "bold"))
ggsave(posname, hm, width = 20, height = 30)