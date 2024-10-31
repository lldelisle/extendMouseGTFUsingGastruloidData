# extend mouse GTF using Gastruloid data

All scripts used to extend the gtf using Gastruloid data before [BREW3R](https://github.com/lldelisle/BREW3R) was polished and released.

The idea is to use both the RefSeq gtf and public RNA-seq full-length stranded and extend the gtf in 3' to be able to get all reads in BRB-seq or 10X sc-RNA-seq.

The repository contains different directories:

- `pipeline` contains all sbatch scripts.
- `scripts` contains R scripts used during the pipeline
- `tables` contains tables with input data etc...
- `outputs` contains the final gtfs that were uploaded to zenodo [here](https://doi.org/10.5281/zenodo.10079673) and [here](https://doi.org/10.5281/zenodo.14016639).

The command lines are written [here](./workflow.sh).


The pipeline is:

- Extend the 'merged_filtered_ensembl' gtf with RefSeq gtf with a custom R script.
- Extend the extended gtf with bulk full length RNA-seq:

  - For each fastq or fastq pair:

    - remove adapters with cutadapt
    - map with STAR with ENCODE parameters
    - get transcripts with Stringtie without reference

  - Combine all Stringtie outputs with Stringtie merge
  - Use the same custom R script to extend the 3' of the input gtf using the output of Stringtie merge.

The new gtf will be extended in 3' end and new exons will be added. This new gtf will loose the specificity of 3'UTR per transcript as it considers gene as the base unit.


Note:

This analysis was performed in January 2023 and unfortunately could not be reproduced in October 2024. One explanation would be that the RefSeq gtf on UCSC website for mm39 changed, I could not find the reason why the results changed for mm10 but it seems to be in the first step during the extension of custom gtf using RefSeq.
