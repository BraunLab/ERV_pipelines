#!/bin/bash
#SBATCH --job-name=star_Index
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

module load STAR/2.7.8a-GCC-10.2.0

STAR --runMode genomeGenerate --runThreadN 4 --limitGenomeGenerateRAM 60000000000  --genomeSAindexNbases 7 --genomeDir . --genomeFastaFiles hervquant_hg19_reference.fa
