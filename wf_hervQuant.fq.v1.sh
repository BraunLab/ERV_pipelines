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

##Path and name for the results folder
out="$path/hervQuant_results"

##Crosscheck whether the folder was created
if [ ! -d "$out" ]; then
	mkdir  -p $out
fi

##Path to folder with all the references
reference=$path/hQ_ref/

##Takes input sample folder name as a command line argument
in_dir=$1

for f in $in_dir/*R1*gz
do
	fR1=$f
done

for f in $in_dir/*R2*gz
do
	fR2=$f
done

##Takes part of the name to name the outputs based on sample names.
out_name=$in_dir

#Makes specific output folder with sample name
mkdir $out/$out_name

##------------------------------------------------------------------------------------------------------------------------------------------

##
# The Log block to print the versions of the program into a log file
##
#############################################################################################################################################

module load STAR/2.7.8a-GCC-10.2.0
module load SAMtools/1.16.1-GCCcore-10.2.0
module load Salmon/1.4.0-gompi-2020b

log="$out/$out_name/$out_name"_program_version.log

echo "Versions log" > $log

echo -e "hervQant piepline version 1\n" >>$log
echo -e "genome-build hg19 masked reference from hervQuant\n" >>$log
echo "STAR Version: 2.7.8a">>$log
echo "Samtools Version: 1.16.1" >>$log
echo "Salmon Version: 1.4.0 " >>$log


##------------------------------------------------------------------------------------------------------------------------------------------
##
# The Workflow block creates from BAM/FASTQ to Salmon Matrix
#############################################################################################################################################

# STAR alignment

STAR \
--genomeDir $reference \
--outFileNamePrefix $out/$out_name/$out_name \
--runThreadN ${SLURM_CPUS_PER_TASK} \
--outFilterMultimapNmax 10 \
--outFilterMismatchNmax 7 \
--readFilesCommand zcat \
--readFilesIn $fR1 $fR2

# filter out all non hERV mapped reads

sam=$out/$out_name/$out_name"Aligned.out.sam" 
bam_filt=$out/$out_name/$out_name".Aligned.sorted.filtered.bam"
sed '/uc.*/d' $sam \
| samtools view -@ ${SLURM_CPUS_PER_TASK} -Shu - \
| samtools sort -@ ${SLURM_CPUS_PER_TASK} -o $bam_filt

# Salmon transcript quantification
time salmon quant \
 -t $reference/hervquant_final_reference.fa \
 -l ISF \
 -a $bam_filt \
 -o $out/$out_name \
 -p ${SLURM_CPUS_PER_TASK}

rm $out/$out_name/$out_name"Aligned.out.sam"

##------------------------------------------------------------------------------------------------------------------------------------------
