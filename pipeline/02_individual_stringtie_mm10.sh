#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 24 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-40 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name Stringtie_indiv # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/new_gtf_mouse/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

# This script run stringtie without reference to get the transcripts from BAM 


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}

### Specify the paths

# Put in dirPathWithResults the directory
# where a directory will be created
# for each sample
dirPathWithResults="$PWD/mm10/"
# All samples are registered into a table where
# first column is the sample name
# We assume that the bam file was generated before
filePathForTable="/home/ldelisle/softwares/extendMouseGTFUsingGastruloidData/tables/table_RNA_all.txt"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve stringtie dependencies
# You can create it with: conda create -n stringtie_R_202211 stringtie r-base  bioconductor-rtracklayer r-plyr
# Comment it if you will use module load
condaEnvName=stringtie_R_202211
######


##################################
####### BEGINING OF SCRIPT #######
##################################

# Check everything is set correctly:
if [ ! -z ${condaEnvName} ]; then
    # Conda environment:
    # This line is to adapt the conda to the shell
    source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
    # We check if the conda environment exists
    exists=$(conda info --envs | awk -v ce=${condaEnvName} '$1==ce{print}' | wc -l)
    # It if does not exists an error is raised
    if [ $exists -ne 1 ]; then
    echo "conda environment ${condaEnvName} does not exists. Create it before."
    exit 1
    fi
    # Activate the conda environment
    conda activate ${condaEnvName}
fi

# Check all softwares are present and write version to stdout:
v=$(stringtie --version)
if [ $? -ne 0 ]
then
  echo "stringtie is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo "stringtie $v"

sample=$(cat ${filePathForTable} | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo $sample

# The analysis part takes part within the pathResults
cd $pathResults

# Check the BAM is present:
if [ ! -e Aligned.sortedByCoord.out.bam ]; then
    echo "The bam file Aligned.sortedByCoord.out.bam is not present run first script first."
    exit 1
fi

if [ ! -e ${sample}_stringtie_min10.gtf ]; then
    stringtie Aligned.sortedByCoord.out.bam \
        -c 10 -s 10 -j 10 \
        -p "${nbOfThreads}" --rf -o ${sample}_stringtie_min10.gtf
fi

mkdir -p ${dirPathWithResults}/allFinalFiles/gtf/

cp *gtf ${dirPathWithResults}/allFinalFiles/gtf/
