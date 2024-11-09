#!/bin/bash

#SBATCH --job-name=New_003                      # Job name
#SBATCH --nodes=1                               # Number of compute nodes
#SBATCH --mem=40G                               # Amount of memory
#SBATCH --partition=standard                    # Partition
#SBATCH --time=60:00:00                         # Wall time
#SBATCH --mail-user=f007qps@dartmouth.edu       # Email 
#SBATCH --mail-type=FAIL                        # Notification type
#SBATCH --output=New_003_%j.out                 # SLURM log file name


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! DEFINE VARIABLES  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WORKINGDIR: path to working directory
# DBNAME: Name of your database directory created in script 001
# EXPERIMEMNTNAME: Prefix for all output files
# REPLICATES: Name of txt file specifying sample and replicate information
# COMPARISONS: Name of txt file specifying DESeq comparison
WORKINGDIR="/dartfs-hpc/rc/home/s/f007qps/final_tRAX_test/tRAX_optimization"
DBNAME="hg38_db"
EXPERIMENTNAME="Opti_Run"
REPLICATES="repfile.txt"
COMPARISONS="design.txt"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! CODE BLOCKS  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- !!!!! Do not edit !!!!! -----#
# Setting paths for code
CODEDIR="$WORKINGDIR/code"
DATABASEDIR="$WORKINGDIR/$DBNAME"
FASTQDIR="$WORKINGDIR/samples"
DATABASEPATH="$DATABASEDIR/db" # Format as <path/to/datavase/db> NOT INTERPRETED AS A FILE OR DIRECTORY.
GTFPATH="$DATABASEDIR/genes.gtf"
REPFILE="$FASTQDIR/$REPLICATES"
PAIRFILE="$FASTQDIR/$COMPARISONS"
BAMDIR="$FASTQDIR/bams/"

# Setting variables for clean up functionality
RESULTSDIR="$FASTQDIR/$EXPERIMENTNAME"
MISMATCHDIR="$FASTQDIR/$EXPERIMENTNAME/mismatch"
THREEPRIME="$MISMATCHDIR/three_prime"
FIVEPRIME="$MISMATCHDIR/five_prime"
DELETIONS="$MISMATCHDIR/deletions"
MISMATCHES="$MISMATCHDIR/mismatches"
COVPLOTS="$FASTQDIR/$EXPERIMENTNAME/coverage_plots"
ENDPOINTS="$FASTQDIR/$EXPERIMENTNAME/endpoint_plots"
DIAGNOSTIC="$FASTQDIR/$EXPERIMENTNAME/Mapping_Diagnostics"
DESEQ="$FASTQDIR/$EXPERIMENTNAME/DESeq_Results"
TAILS="$FASTQDIR/$EXPERIMENTNAME/tRNA_tails"
FINALOUT="$WORKINGDIR/Outputs"


# ----- Start
date=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job started at: $date"

# ----- Move into fastq directory
if [[ -d "$FASTQDIR" && "$(ls -A "$FASTQDIR")" ]]; then
    echo "Sample directory exists and is not empty. Proceeding..."
    cd "$FASTQDIR"
else
    echo "Sample directory either does not exist or is empty. Exiting."
    exit 1
fi

# ----- Check that the replicate file exists
if [ ! -e "$REPFILE" ]; then
    echo "Cannot find '$REPFILE'. Exiting"
    exit 3
else 
    echo "'$REPFILE' exists. Continuing..."
fi

# ----- Check that the pair file exists
if [ ! -e "$PAIRFILE" ]; then
    echo "Cannot find '$PAIRFILE'. Exiting"
    exit 3
else 
    echo "'$PAIRFILE' exists. Continuing..."
fi

# ----- Check that the gtf file exists
if [ ! -e "$GTFPATH" ]; then
    echo "Cannot find '$GTFPATH'. Exiting"
    exit 3
else 
    echo "'$GTFPATH' exists. Continuing..."
fi

# ----- Output directory generation
if [ ! -d "$BAMDIR" ]; then
    echo "Creating BAM directory..."
    mkdir -p "$BAMDIR"
    echo "BAM directory successfully created."
else
    echo "Directory already exists. Continuing..."
fi

# ----- Activate conda environment
echo "Beginning mapping and downstream analysis"
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env

# ----- Run 03_processsamples.py 
python "$CODEDIR"/03_processsamples.py \
--experimentname="$EXPERIMENTNAME" \
--databasename="$DATABASEPATH" \
--ensemblgtf="$GTFPATH" \
--samplefile="$REPFILE" \
--exppairs="$PAIRFILE" \
--bamdir="$BAMDIR" \
--nofrag 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! DIRECTORY CLEANING  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ----- Move into mismatch directory
if [[ -d "$MISMATCHDIR" && "$(ls -A "$MISMATCHDIR")" ]]; then
    echo "mismatch exists and is not empty. Proceeding..."
    cd "$MISMATCHDIR"
else
    echo "mismatch directory either does not exist or is empty. Exiting."
    exit 1
fi

#----- Remove all files related to individual positions of tRNA reads
PREDUMP=$(ls | wc -l)
echo "There are $PREDUMP files in the directory '$MISMATCHDIR'."
echo "Removing unnecessary files"
rm *_posaminodeletions.pdf
rm *_posaminomismatches.pdf
rm *_possamplereadstarts.pdf
rm *_possampledeletions.pdf
rm *_possamplemismatches.pdf
POSTDUMP=$(ls | wc -l)
REMOVED=$((PREDUMP-POSTDUMP))
echo "$REMOVED files were removed."

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Output directory generation
if [ ! -d "$THREEPRIME" ]; then
    echo "Creating 3' mismatch directory..."
    mkdir -p "$THREEPRIME"
    echo "3' directory successfully created."
else
    echo "3' directory already exists. Continuing..."
fi

# Move all 3' related files into the 3' directory
mv *_threeprimecounts.pdf "$THREEPRIME"
mv *_threeprimeheatmap.pdf "$THREEPRIME"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Output directory generation
if [ ! -d "$FIVEPRIME" ]; then
    echo "Creating 5' mismatch directory..."
    mkdir -p "$FIVEPRIME"
    echo "5' directory successfully created."
else
    echo "5' directory already exists. Continuing..."
fi

# Move all 5' related files into the 3' directory
mv *_fiveprimecounts.pdf "$FIVEPRIME"
mv *_fiveprimeheatmap.pdf "$FIVEPRIME"
mv *_fiveprimeends.pdf "$FIVEPRIME"
mv *-trnapositionfiveprime.pdf "$FIVEPRIME"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Output directory generation
if [ ! -d "$DELETIONS" ]; then
    echo "Creating deletions mismatch directory..."
    mkdir -p "$DELETIONS"
    echo "Deletions directory successfully created."
else
    echo "Deletions directory already exists. Continuing..."
fi

# ----- Move all deletion related files into the 3' directory
mv *_delete.pdf "$DELETIONS"
mv *_deletionheatmap.pdf "$DELETIONS"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Output directory generation
if [ ! -d "$MISMATCHES" ]; then
    echo "Creating mismatches mismatch directory..."
    mkdir -p "$MISMATCHES"
    echo "Mismatches directory successfully created."
else
    echo "Mismatches directory already exists. Continuing..."
fi

# ----- Move files into results directory
mv *_mismatch.pdf "$MISMATCHES"
mv *_mismatchheatmap.pdf "$MISMATCHES"
mv *-trnapositionmismatches.pdf "$MISMATCHES"
mv *-Full_relative_summed_tRNA_mismatches.png "$RESULTSDIR"

# ----- Move into root output directory
if [[ -d "$RESULTSDIR" && "$(ls -A "$RESULTSDIR")" ]]; then
    echo "Results directory exists and is not empty. Proceeding..."
    cd "$RESULTSDIR"
else
    echo "Results directory either does not exist or is empty. Exiting."
    exit 1
fi

# ----- Output directory generation
if [ ! -d "$COVPLOTS" ]; then
    echo "Creating coverage plot directory..."
    mkdir -p "$COVPLOTS"
    echo "Coverage plot directory successfully created."
else
    echo "Coverage plot directory already exists. Continuing..."
fi

# ----- Move all coverage plots into the coverage plot directory
echo "Moving all coverage plots to '$COVPLOTS'."
mv *_cov.pdf "$COVPLOTS"
mv *-combinedcoverages.pdf "$COVPLOTS"
mv *-coverage.pdf "$COVPLOTS"

# ----- Output directory generation
if [ ! -d "$ENDPOINTS" ]; then
    echo "Creating endpoint plot directory..."
    mkdir -p "$ENDPOINTS"
    echo "Endpoint plot directory successfully created."
else
    echo "Endpoint plot directory already exists. Continuing..."
fi

# ----- Move all coverage plots into the coverage plot directory
echo "Moving all endpoint plots to '$ENDPOINTS'."
mv *_endplot.pdf "$ENDPOINTS"
mv *_trnaendcounts.txt "$ENDPOINTS"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Output directory generation
if [ ! -d "$DIAGNOSTIC" ]; then
    echo "Creating diagnostic information directory..."
    mkdir -p "$DIAGNOSTIC"
    echo "Diagnostic information directory successfully created."
else
    echo "Diagnostic information directory already exists. Continuing..."
fi

mv *-mapinfo.pdf "$DIAGNOSTIC"
mv *-trnamapinfo.pdf "$DIAGNOSTIC"
mv *-trnamapinfo.txt "$DIAGNOSTIC"
mv *-mapinfo.txt "$DIAGNOSTIC"
mv *-mapstats.txt "$DIAGNOSTIC"
mv *-runinfo.txt "$DIAGNOSTIC"
mv *-qa.html "$DIAGNOSTIC"
mv *-aminocounts.pdf "$DIAGNOSTIC"
mv *-aminorealcounts.pdf "$DIAGNOSTIC"
mv *-coverages.txt "$DIAGNOSTIC"
mv *-coverages.pdf "$DIAGNOSTIC"
mv *-readlengths.pdf "$DIAGNOSTIC"
mv *-readlengths.txt "$DIAGNOSTIC"
mv *-typerealcounts.pdf "$DIAGNOSTIC"
mv *-typerealcounts.txt "$DIAGNOSTIC"
mv *-mismatches.txt "$DIAGNOSTIC"
mv *-mismatches.pdf "$DIAGNOSTIC"
mv *-typecounts.pdf "$DIAGNOSTIC"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ----- Output directory generation
if [ ! -d "$DESEQ" ]; then
    echo "Creating DESeq results directory..."
    mkdir -p "$DESEQ"
    echo "DESeq results directory successfully created."
else
    echo "DESeq results directory already exists. Continuing..."
fi

mv *-aminoscatter.pdf "$DESEQ"
mv *-countcompare.pdf "$DESEQ"
mv *-typescatter.pdf "$DESEQ"
mv *-volcano_tRNA.pdf "$DESEQ"
mv *-volcano.pdf "$DESEQ"
mv *-pca.pdf "$DESEQ"
mv *-pcatrna.pdf "$DESEQ"
mv *-lengthcompare.pdf "$DESEQ"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if [ ! -d "$DESEQ/tabular_data" ]; then
    echo "Creating DESeq tabular data subdirectory..."
    mkdir -p "$DESEQ/tabular_data"
    echo "DESeq tabular data successfully created."
else
    echo "DESeq tabular data already exists. Continuing..."
fi

mv *-dispersions.txt $DESEQ/tabular_data
mv *-SizeFactors.txt $DESEQ/tabular_data
mv *-logvals.txt $DESEQ/tabular_data
mv *-medians.txt $DESEQ/tabular_data
mv *-normalizedreadcounts.txt $DESEQ/tabular_data
mv *-padjs.txt $DESEQ/tabular_data
mv *-readcounts.txt $DESEQ/tabular_data
mv *-combine.txt $DESEQ/tabular_data
mv *-trnacounts.txt $DESEQ/tabular_data
mv *-aminocounts.txt $DESEQ/tabular_data
mv *-anticodoncounts.txt $DESEQ/tabular_data
mv *-genetypes.txt $DESEQ/tabular_data
mv *-typecounts.txt $DESEQ/tabular_data
mv *-typecountsreal.txt $DESEQ/tabular_data

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cd "$RESULTSDIR"
if [ ! -d "$TAILS" ]; then
    echo "Creating tRNA_tails directory..."
    mkdir -p "$TAILS"
    echo "tRNA_tails directory successfully created."
else
    echo "tRNA_tails directory already exists. Continuing..."
fi

mv *-trnaendcounts.txt "$TAILS"
mv *-trnatails.pdf "$TAILS"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ----- Back up one directory to the samples
cd "$FASTQDIR"

# ----- Move final files into resultsDir
mv mismatchcompare.txt "$RESULTSDIR"
mv positiondeletions.txt "$RESULTSDIR"
mv positionmismatches.txt "$RESULTSDIR"
mv Rplots.pdf positionmismatches.pdf 
mv positionmismatches.pdf "$MISMATCHES"

# ----- Move all output directories to finalOut
cd $WORKINGDIR
if [ ! -d "$FINALOUT" ]; then
    echo "Creating final output directory..."
    mkdir -p "$FINALOUT"
    echo "Final output directory successfully created."
else
    echo "Final output directory already exists. Continuing..."
fi

# ----- Moving straggling files to relevant outputs
mv "$RESULTSDIR/mismatchcompare.txt" "$MISMATCHES"
mv "$RESULTSDIR/positionmismatches.pdf" "$MISMATCHES"
mv "$RESULTSDIR/positionmismatches.txt" "$MISMATCHES"
mv "$RESULTSDIR/positiondeletions" "$DELETIONS"
mv "$RESULTSDIR" "$FINALOUT"

# -----End
endDate=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job ended at: $endDate"
