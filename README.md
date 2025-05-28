# ERV\_pipelines

## Overview

This repository provides scripts, reference files, and conda environments required to quantify endogenous retrovirus (ERV) expression from RNA-Seq data. It supports analysis from either BAM or FASTQ input files.

The repository includes two distinct ERV quantification workflows:

* **hQ (hervQuant)**: From the Vincent Lab – [hervQuant GitHub](https://unclineberger.org/vincentlab/resources/)
* **eNM (ERV NatMed)**: Based on the methodology described in [Braun et al., *Nature Medicine*, 2020](https://www.nature.com/articles/s41591-020-0839-y)

---

## Directory Structure

```
ERV_pipelines/
├── eNM_ref/                  # Reference files for eNM pipeline (hg19, Bowtie2 + HTSeq)
├── erv-pipe.yml             # Conda environment YAML for eNM pipeline
├── hQ_ref/                  # Reference files for hQ pipeline (hg19, STAR + Salmon)
├── README.md
├── wf_erv_natmed.bam.v1.sh  # eNM pipeline (BAM input)
├── wf_erv_natmed.fq.v1.sh   # eNM pipeline (FASTQ input)
├── wf_hervQuant.bam.v1.sh   # hQ pipeline (BAM input)
├── wf_hervQuant.fq.v1.sh    # hQ pipeline (FASTQ input)
└── wf_joblist.v1.sh         # Submission helper script
```

---

## Reference Directories

### `hQ_ref/`

Contains all reference files based on **hg19**, required to:

* Align reads using **STAR**
* Quantify expression using **Salmon**

### `eNM_ref/`

Contains all reference files based on **hg19**, required to:

* Align reads using **Bowtie2**
* Quantify expression using **HTSeq**

---

## Conda Environment Setup

To run the **eNM** pipeline, set up the conda environment as follows:

```bash
conda env create -f erv-pipe.yml --name erv-pipe
```

Activate it before running the pipeline:

```bash
conda activate erv-pipe
```

---

## Workflow Scripts

Choose the appropriate script depending on whether your input is in **BAM** or **FASTQ** format.

### eNM Pipeline

* `wf_erv_natmed.bam.v1.sh`
* `wf_erv_natmed.fq.v1.sh`

### hQ Pipeline

* `wf_hervQuant.bam.v1.sh`
* `wf_hervQuant.fq.v1.sh`

---

## Running the Pipelines

### 1. Set Up

Download or clone the repository into the same folder containing your BAM or FASTQ files.

Repository location on Farnam (Yale HPC):
`/gpfs/gibbs/pi/braun/workflows/ERV`

GitHub:
[https://github.com/BraunLab/ERV\_pipelines](https://github.com/BraunLab/ERV_pipelines)

---

### 2. Submit the Job List

```bash
sbatch wf_joblist.v1.sh
```

This generates a submission script (e.g., `dsq-joblist-2022-09-14.sh`) using [dSQ](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/) for batch scheduling.

You will see a message indicating the generated submission file. Then run:

```bash
sbatch dsq-joblist-2022-09-14.sh
```

---

### 3. Monitor Job Progress

```bash
module load dSQ
dsqa -j <jobid>
```

Replace `<jobid>` with your actual job ID.

---

### 4. Rerun Failed Jobs (if any)

If jobs fail due to resource issues:

```bash
module load dSQ
dsqa -j <jobid> -f jobsfile.txt > re-run_jobs.txt 2> <jobid>_report.txt

dsq --job-file re-run_jobs.txt --time=2-00:00:00 --cpus-per-task=5 --mem=35G --mail-type=ALL
```

For more information, visit the [Yale dSQ documentation](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/).

---

## Questions or Support

For questions, feedback, or assistance, please contact:
📧 **deepak(dot)poduval(at)yale(dot)edu**
