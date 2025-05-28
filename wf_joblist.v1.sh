#!/bin/bash

## Loop through bam files / fastq file containing folders and create joblist to make dSQ submission. 
## The script produces 2 dsq joblist that can be submitted as sbatch dsq-joblist-yyyy-mm-dd.sh

# Replace the prefix of bam accordingly to the project. For fastq file provide the prefix for fastq folder names
# Replace use appropriate hervquant and erv natmed pipeline 

for file in ACH-*bam
do
	echo "bash wf_hervQuant.bam.v1.sh $file" >> joblist_hQ.txt
done

for file in ACH-*bam
do
	echo "bash wf_erv_natmed.bam.v1.sh $file" >> joblist_eNM.txt
done

## Load dSQ module and create serial job submission
module load dSQ
dsq --job-file joblist_hQ.txt -p day --time=12:00:00 --cpus-per-task=5 --mem=35G --mail-type=ALL
dsq --job-file joblist_eNM.txt -p day --time=12:00:00 --cpus-per-task=5 --mem=35G --mail-type=ALL

#--------------------------------------------------------------------------------------------
