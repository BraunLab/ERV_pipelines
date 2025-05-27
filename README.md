# ERV_pipelines

## README

## TL;DR

This directory includes a GTEx (Release V10) singularity container of scripts and executables, and 2 folders of GTEx pipeline to be run on either bam or fastq file named accordingly. Each folder contains a workflow bash script and two slurm scripts. This need to be run inside the directory containing experiment files. The bam pipeline using GRCh38 human and bam files is tested and completed on the Farnam HPC cluster at Yale.

## Directory

### references

The directory contains a genome fasta file, indices, and genes gtf file for human version GRCh38.p10 (NCBI:GCA_000001405.25). A gene model for the gtf file was made using a script from GTEx. Two other folders of reference indices one for STAR, mapping reads, and RSEM, quantification. This gives the flexibility to change references to a newer version or a different version for other purposes.

## Singularity container

### gtex_rnaseq_v10.sif

This is the singularity container that includes all the necessary helper scripts and programs that need to run the RNA-Seq analysis including STAR, Picard, RNA-SeQC2, and RSEM. The container also has helper scripts to combine the results from individual runs.

**Note: Copy the references and singularity container from "/SAY/standard/braunlab-CC1413-MEDCCC/BraunLab_pipelines/gtex_rnaseq_pipeline_v1". Also, the appropriate folder of scripts is present at the same locations. A copy of the scripts is also available at <https://github.com/BraunLab/gtex_rnaseq_pipeline_v1>.**

## Bash script

### wf_gtex_rnaseq.v1.sh

The bash script contains sequential steps required to do the RNA-Seq analysis from bam to count matrix files for individual bam files. This script takes the input bam file and runs all essential steps including creating an output log file, for the genome version and program version that was run.

## Slurm scripts

Slurm is a scheduler used in Farnam to schedule jobs.

### wf_gtex_rnaseq_joblist.v1.slurm

Farnam uses a dSQ scheduler to serialize jobs. This slurm script will loop through the names of bam files making a joblist txt file, with each line containing the execution of workflow bash script for each bam file. It will submit the joblist to the dSQ scheduler outputting the submission line to run the whole pipeline in serial.

### wf_gtex_rnaseq_aggregate.v1.slurm

This slurm script is used to make result aggregate from RSEM and RNASeQC2 runs using helper scripts.

## Instructions

Download this directory's contents to your destination folder with all the bam files to analyze.
The first step in the workflow is to run

```bash
sbatch wf_gtex_rnaseq_joblist.v1.slurm
```

This step will create a joblist file and submit it to dSQ serial job submission script with a date of the run.
e.g. dsq-joblist-2022-09-14.sh
It will output a message for the exact command to run.

```bash
sbatch dsq-joblist-2022-09-14.sh
```

This block of code will run the workflow on all the bam files with available resources in Farnam with a "jobid". To check the progress of the run the code block is below. It will give you the status of the progress.

```bash
module load dSQ
dsqa -j <jobid>
```

If the run fails due to resource unavailability. If you haven't loaded the dSQ module load it again and run the code block below. That should rerun all the jobs that failed.
Please refer to <https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/> to learn more about dSQ.

```bash
#module load dSQ
dsqa -j <jobid> -f jobsfile.txt > re-run_jobs.txt 2> <jobid>_report.txt
dsq --job-file re-run_jobs.txt --time=2-00:00:00 --cpus-per-task=5 --mem=35G --mail-type=ALL
```

Please feel free to contact me for any further queries  **deepak(dot)poduval(at)yale(dot)edu**
