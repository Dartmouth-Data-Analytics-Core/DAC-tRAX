#!/bin/bash

#SBATCH --job-name=New_001                         # Job name (!!! EDIT ME !!!)
#SBATCH --nodes=1                                  # Number of compute nodes
#SBATCH --mem=40G                                  # Amount of memory
#SBATCH --partition=standard                       # Partition
#SBATCH --time=60:00:00                            # Wall time
#SBATCH --mail-user=f007qps@dartmouth.edu          # Email (!!! EDIT ME !!!)
#SBATCH --mail-type=FAIL                           # Notification type
#SBATCH --output=New_001_%j.out                    # SLURM log file name (!!! EDIT ME !!!)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! DEFINE VARIBALES  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WORKINGDIR: Absolute path to your working directory
# GENOME: One of the following options: mm10, mm10mito, rn6, hg19, hg19mito, hg38, hg38mito, sacCer3, dm6
WORKINGDIR="/dartfs-hpc/rc/home/s/f007qps/final_tRAX_test/tRAX_optimization"
GENOME="hg38"





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~! CODE BLOCKS  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----- !!!!! Do not edit !!!!! -----#
CODEDIR="$WORKINGDIR/code"
DATABASEDIR="$WORKINGDIR/${GENOME}_db"

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
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/trax_env

# ----- Check if conda activation was successful
if [[ $? -ne 0 ]]; then
    echo "Failed to activate conda environment. Exiting."
    exit 3
else
    echo "trax_env successfully activated."
fi

# ----- Output directory generation
if [ ! -d "$DATABASEDIR" ]; then
    echo "Creating database directory..."
    mkdir -p "$DATABASEDIR"
    echo "Database directory successfully created."
else
    echo "Directory already exists. Continuing..."
fi

# ----- Run code to generate database
cd "$CODEDIR"
echo "Starting 01_quickdb.bash..."
bash 01_quickdb.bash "$GENOME" "$DATABASEDIR"

# -----End
endDate=$(date '+%Y-%m-%d %H:%M:%S')
echo "Job ended at: $endDate"


