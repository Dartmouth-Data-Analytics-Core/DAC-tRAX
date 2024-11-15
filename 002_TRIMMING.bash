#!/bin/bash

#SBATCH --job-name=New_002                             # Job name
#SBATCH --nodes=1                                  # Number of compute nodes
#SBATCH --mem=40G                                  # Amount of memory
#SBATCH --partition=standard                       # Partition
#SBATCH --time=60:00:00                            # Wall time
#SBATCH --mail-user=f007qps@dartmouth.edu          # Email 
#SBATCH --mail-type=FAIL                           # Notification type
#SBATCH --output=New_002_%j.out                        # SLURM log file name

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! DEFINE VARIBALES  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WORKINGDIR: Absolute path to working directory
# RUNFILE: Name of your runfile.txt file
# LAYOUT: Either SE for single-end and PE for paired-end.
WORKINGDIR="/dartfs-hpc/rc/home/s/f007qps/final_tRAX_test/tRAX_optimization"
RUNFILE="runfile.txt"
LAYOUT="SE"




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! CODE BLOCKS  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- !!!!! Do not edit !!!!! -----#
CODEDIR="$WORKINGDIR/code"
FASTQDIR="$WORKINGDIR/samples"
OPNAME="trim"
SAMPLEFILE="$FASTQDIR/$RUNFILE"
DIAGNOSTICS="$FASTQDIR/Trimming_Diagnostics"

# ----- Start
date=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job started at: $date"

# ----- Move into the directory that holds all the tRAX code
echo "Moving to working directory..."

# ----- Check that the working directory exists and is not empty
if [[ -d "$WORKINGDIR" && "$(ls -A "$WORKINGDIR")" ]]; then
    echo "Working directory exists and is not empty. Proceeding..."
    cd "$WORKINGDIR"
else
    echo "Working directory either does not exist or is empty. Exiting."
    exit 1
fi

# ----- Load Conda
source /optnfs/common/miniconda3/etc/profile.d/conda.sh

# ----- Check that conda was sucessfully sources
if ! command -v conda &> /dev/null; then
    echo "Failed to source conda. Exiting."
    exit 2
else    
    echo "Successfully sourced conda"
fi

# ----- Activate conda environment
echo "Activing conda environment: trax_env"
conda activate /dartfs/rc/lab/O/OrellanaE/GMBSR/envs/trax_env

# ----- Check if conda activation was successful
if [[ $? -ne 0 ]]; then
    echo "Failed to activate conda environment. Exiting."
    exit 3
else
    echo "trax_env successfully activated."
fi

# ----- Move into sample directory
if [[ -d "$FASTQDIR" && "$(ls -A "$FASTQDIR")" ]]; then
    echo "Sample directory exists and is not empty. Proceeding..."
    cd "$FASTQDIR"
else
    echo "Sample directory either does not exist or is empty. Exiting."
    exit 1
fi

# ----- Check that the runfile exists
if [ ! -e "$SAMPLEFILE" ]; then
    echo "Cannot find '$SAMPLEFILE'. Exiting"
    exit 4
else 
    echo "'$SAMPLEFILE' exists. Continuing..."
fi

# ----- Run code to trim adapters in either SE or PE mode
if [ "$LAYOUT" == "SE" ]; then
    echo "Beginning adapter trimming in SE mode"
    python "$CODEDIR"/02_trimadapters.py \
    --runname "$OPNAME" \
    --runfile "$SAMPLEFILE" \
    --singleend
else
    echo "Beginning adapter trimming in PE mode"
    python "$CODEDIR"/02_trimadapters.py \
    --runname "$OPNAME" \
    --runfile "$SAMPLEFILE"
fi
echo "Adapter trimming completed!"

# ----- Move into the sample directory
cd "$FASTQDIR"

# Create trimming diagnostics
if [ ! -d "$DIAGNOSTICS" ]; then
    echo "Creating Trim diagnostics directory..."
    mkdir -p "$DIAGNOSTICS"
    echo "Trim diagnostics directory successfully created."
else
    echo "Trim diagnostics directory already exists. Continuing..."
fi

# ----- Move diagnostic files
mv trim_ca.pdf "$DIAGNOSTICS"
mv trim_ca.txt "$DIAGNOSTICS"
mv trimindex.txt "$DIAGNOSTICS"
mv trim_log.txt "$DIAGNOSTICS"
mv trim_manifest.txt "$DIAGNOSTICS"

# ----- Rename the Rplots and move
mv Rplots.pdf Trimming_percentges.pdf
mv Trimming_percentges.pdf "$DIAGNOSTICS"

# -----End
endDate=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job ended at: $endDate"


