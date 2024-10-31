#!/bin/bash

# Specify location of this directory:
gitHubDirectory=/home/ldelisle/softwares/extendMouseGTFUsingGastruloidData/

# This has been done on helvetios
# The first thing to do is to create a working directory:
wd=/scratch/ldelisle/new_gtf_mouse/
mkdir -p ${wd}
cd ${wd}

# Get the index from work on scratch:
cp -r /work/updub/scratch/ldelisle/genomes/STARIndex_2.7.9a/mm10/ /scratch/ldelisle/genomes/STARIndex_2.7.9a/

mkdir -p ${wd}/mm10/
# Get the gtf from mm10 from zenodo:
wget https://zenodo.org/record/7510406/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz?download=1 -O ${wd}/mm10/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz
gunzip ${wd}/mm10/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf.gz

# Get the refseq:
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
mv mm10.ncbiRefSeq.gtf.gz ${wd}/mm10/
gunzip -k ${wd}/mm10/mm10.ncbiRefSeq.gtf.gz

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# conda create -n stringtie_R_202211 stringtie r-base bioconductor-rtracklayer
# r-base version 4.2.2
# rtracklayer version 1.58.0
# stringtie version 2.2.1
conda activate stringtie_R_202211


# Extend the gtf with Refseq:
Rscript ${gitHubDirectory}/scripts/improve_gtf3p.R ${wd}/mm10/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC.gtf \
    ${wd}/mm10/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsCDSOnly_UCSC_Refseq.gtf \
    ${wd}/mm10/mm10.ncbiRefSeq.gtf.gz

conda deactivate

# Run the SR:
sbatch ${gitHubDirectory}/pipeline/01_RNAseq_SR_mm10.sh
# Run the PE:
sbatch ${gitHubDirectory}/pipeline/01_RNAseq_PE_mm10.sh

filePathForTable="/home/ldelisle/softwares/extendMouseGTFUsingGastruloidData/tables/table_RNA_all.txt"

cat ${gitHubDirectory}/tables/table_RNA_SR.txt ${gitHubDirectory}/tables/table_RNA_PE.txt | cut -f 1 > ${filePathForTable}

sbatch ${gitHubDirectory}/pipeline/02_individual_stringtie_mm10.sh

sbatch ${gitHubDirectory}/pipeline/03_extend_gtf_mm10.sh

gzip -c ${wd}/mm10/allFinalFiles/gtf/mm10_custom102_allGastruloids_min10_extended.gtf > ${gitHubDirectory}/outputs/mm10_custom102_allGastruloids_min10_extended.gtf.gz
cp ${wd}/mm10/mm10_custom102_allGastruloids_min10/mm10_custom102_allGastruloids_min10_*log ${gitHubDirectory}/outputs/

##################################################
# mm39
##################################################

# Get the index from work on scratch:
cp -r /work/updub/scratch/ldelisle/genomes/STARIndex_2.7.9a/mm39/ /scratch/ldelisle/genomes/STARIndex_2.7.9a/

mkdir -p ${wd}/mm39/
# Get the gtf from mm39 from zenodo:
wget https://zenodo.org/record/7510797/files/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf.gz?download=1 -O ${wd}/mm39/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf.gz
gunzip ${wd}/mm39/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf.gz

# Get the refseq:
# The file currently on the webserver is different from the file downloaded in 2023-10
# It was uploaded in 2023-08-29
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz
mv mm39.ncbiRefSeq.gtf.gz ${wd}/mm39/
gunzip -k ${wd}/mm39/mm39.ncbiRefSeq.gtf.gz

# Extend the gtf with Refseq:
conda activate stringtie_R_202211
Rscript ${gitHubDirectory}/scripts/improve_gtf3p.R ${wd}/mm39/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC.gtf \
    ${wd}/mm39/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.108_ExonsCDSOnly_UCSC_Refseq.gtf \
    ${wd}/mm39/mm39.ncbiRefSeq.gtf.gz
# Found 45081 potential exon extensions. Will now check for collision.
# Found 2760 exon extensions which collide.
# Among them 2698 exon extensions overlap a larger region.
# Therefore only 43509 exons will be extended.
# Found 969 transcripts to extend.
# Added 1564 exons.

conda deactivate

# Create symlink for cutadapt results:
for d in  mm10/Rekaik_wt_* mm10/Beccari*; do
    sample=$(basename $d)
    mkdir -p mm39/$sample
    ln -s $PWD/${d}/*cutadapt* mm39/$sample/
done
# Run the SR and PE:
sbatch ${gitHubDirectory}/pipeline/01_RNAseq_SR_mm39.sh
sbatch ${gitHubDirectory}/pipeline/01_RNAseq_PE_mm39.sh

sbatch --dependency=afterok:27838625,27838626 ${gitHubDirectory}/pipeline/02_individual_stringtie_mm39.sh
sbatch --dependency=afterok:27839181 ${gitHubDirectory}/pipeline/03_extend_gtf_mm39.sh

gzip -c ${wd}/mm39/allFinalFiles/gtf/mm39_custom108_allGastruloids_min10_extended.gtf > ${gitHubDirectory}/outputs/mm39_custom108_allGastruloids_min10_extended.gtf.gz
cp ${wd}/mm39/mm39_custom108_allGastruloids_min10/mm39_custom108_allGastruloids_min10_*log ${gitHubDirectory}/outputs/
