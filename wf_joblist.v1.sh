#!/bin/bash

#############################################################################################
## Loop through bam files and create joblist to make dSQ submission

for file in ACH-*bam
do
	echo "bash wf_erv_pipeline.sh $file" >> joblist_erv_natmed.txt
done

## Load dSQ module and create serial job submission
module load dSQ
dsq --job-file joblist_erv_natmed.txt -p day --time=12:00:00 --cpus-per-task=5 --mem=35G --mail-type=ALL

#--------------------------------------------------------------------------------------------
