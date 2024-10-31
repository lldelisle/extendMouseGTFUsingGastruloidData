#!/bin/bash

#SBATCH -o slurm-%x-%A.out # Template for the std output of the job uses the job name and the job id
#SBATCH -e slurm-%x-%A.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 24 # This allows to speed the mapping part of the pipeline
#SBATCH --time 04:00:00 # This depends on the size of the fastqs
#SBATCH --job-name Stringtie_merge # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/new_gtf_mouse/ # This directory must exist, this is where will be the error and out files
## Specific to baobab:
##SBATCH --partition=shared-cpu # shared-cpu for CPU jobs that need to run up to 12h, public-cpu for CPU jobs that need to run between 12h and 4 days

# This script run stringtie merge without reference


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
dirPathWithResults="$PWD/mm39/"
# The Rscript improve_gtf3p.R must be in
dirPathForScripts="/home/ldelisle/softwares/extendMouseGTFUsingGastruloidData/scripts/"
filePathForGTF="${dirPathWithResults}/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC_Refseq.gtf"

### Specify the way to deal with dependencies:

# The module load depends on the installation
# # For baobab:
# module purge
# module load Miniconda3/4.9.2

# You can choose to use a conda environment to solve stringtie dependencies
# You can create it with: conda create -n stringtie_R_202211 stringtie r-base bioconductor-rtracklayer r-plyr
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
    # source /home/ldelisle/miniconda2/etc/profile.d/conda.sh

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
R --version
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
# Check rtracklayer and plyr are installed:
Rscript -e "library(rtracklayer);library(plyr)"
if [ $? -ne 0 ]
then
  echo "Some R packages are missing check rtracklayer and plyr are installed."
  exit 1
fi
# Check the Rscript:
if [ ! -e ${dirPathForScripts}/improve_gtf3p.R ]; then
    echo "The Rscript improve_gtf3p.R is not in ${dirPathForScripts}/"
    exit 1
fi

sample="mm39_custom108_allGastruloids_min10"

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo $sample

# The analysis part takes part within the pathResults
cd $pathResults

if [ ! -e ${sample}_stringtie_merge.gtf ]; then
    stringtie --merge -p ${nbOfThreads} \
    -o ${sample}_stringtie_merge.gtf \
    ${dirPathWithResults}/allFinalFiles/gtf/*_stringtie_min10.gtf
fi


if [ ! -e ${sample}_extended.gtf ]; then
    Rscript ${dirPathForScripts}/improve_gtf3p.R \
    $filePathForGTF \
    ${sample}_extended.gtf \
    ${sample}_stringtie_merge.gtf > ${sample}_extension.log
fi

wget "https://raw.githubusercontent.com/lldelisle/toolBoxForMutantAndWTGenomes/main/scripts/convertGtfToRefGeneForIGV.R" \
  -O ${dirPathForScripts}/convertGtfToRefGeneForIGV.R -nc

if [ ! -e ${sample}_RefGene.txt ]; then
  Rscript ${dirPathForScripts}/convertGtfToRefGeneForIGV.R ${sample}_extended.gtf
  mv RefGene.txt ${sample}_RefGene.txt
fi

mkdir -p ${dirPathWithResults}/allFinalFiles/gtf/

cp *gtf ${dirPathWithResults}/allFinalFiles/gtf/

cp *_RefGene.txt ${dirPathWithResults}/allFinalFiles/gtf/
