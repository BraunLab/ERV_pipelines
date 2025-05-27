#!/bin/bash -ue

# wf_hervQuant.v1.sh 

###########################################################################################################################################
# In-house hervQuant pipeline to analyse RNA-Seq: Version 1										  #
#																	  #
###########################################################################################################################################

##
# The Path block to setup path of folders and files - input and outputs
############################################################################################################################################

##Path of the working folder
path=`pwd`

in=$1
n=${in%.A*}

##Path and name for the results folder
out="$path/hervQuant_results"

##Crosscheck whether the folder was created
if [ ! -d "$out" ]; then
	mkdir  -p $out
fi

mkdir $out/$n

##Path to folder with all the references
reference=$path/hQ_ref/
##------------------------------------------------------------------------------------------------------------------------------------------

##
# The Log block to print the versions of the program into a log file
##
#############################################################################################################################################

module load GATK/4.4.0.0-GCCcore-10.2.0-Java-17
module load STAR/2.7.8a-GCC-10.2.0
module load SAMtools/1.16.1-GCCcore-10.2.0
module load Salmon/1.4.0-gompi-2020b
module load pigz/2.6-GCCcore-10.2.0

log=$out/$n/$n"_program_version.log"

echo "Versions log" > $log

echo -e "hervQant piepline version 1\n" >>$log
echo -e "genome-build hg19 masked reference from hervQuant\n" >>$log
echo "GATK0 Version: 2.7.8a">>$log
echo "STAR Version: 2.7.8a">>$log
echo "Samtools Version: 1.16.1" >>$log
echo "Salmon Version: 1.4.0 " >>$log

##------------------------------------------------------------------------------------------------------------------------------------------
##
# The Workflow block creates from BAM/FASTQ to Salmon Matrix
#############################################################################################################################################

time samtools index $in

time gatk --java-options "-Xmx35G -XX:ParallelGCThreads=5" SamToFastq \
       -I $in \
       -NON_PF true \
       --VALIDATION_STRINGENCY SILENT \
       -F $out/$n/$n"_1.fastq" \
       -F2 $out/$n/$n"_2.fastq" 

time pigz $out/$n/$n*".fastq"

fR1=$out/$n/$n"_1.fastq.gz"
fR2=$out/$n/$n"_2.fastq.gz"

# STAR alignment

time STAR \
--genomeDir $reference \
--outFileNamePrefix $out/$n/$n \
--runThreadN ${SLURM_CPUS_PER_TASK} \
--outFilterMultimapNmax 10 \
--outFilterMismatchNmax 7 \
--readFilesCommand zcat \
--readFilesIn $fR1 $fR2

rm $fR1 $fR2

# filter out all non hERV mapped reads

sam=$out/$n/$n"Aligned.out.sam" 
bam_filt=$out/$n/$n".Aligned.sorted.filtered.bam"
sed '/uc.*/d' $sam \
| samtools view -@ ${SLURM_CPUS_PER_TASK} -Shu - \
| samtools sort -@ ${SLURM_CPUS_PER_TASK} -o $bam_filt

rm $sam

# Salmon transcript quantification
time salmon quant \
 -t $reference/hervquant_final_reference.fa \
 -l ISF \
 -a $bam_filt \
 -o $out/$n \
 -p ${SLURM_CPUS_PER_TASK}

##------------------------------------------------------------------------------------------------------------------------------------------
