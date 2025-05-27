#!/bin/bash -ue

module purge

ml miniconda Java/1.8.345

conda activate erv-pipe

##Path of the working folder
#path="/home/dp935/palmer_scratch/ERV/DepMap_CCLE"
path=`pwd`
#bam=$1

# Check if the folder name is provided as an argument
if [ $# -eq 0 ]; then
    echo "No folder name provided. Please provide the folder name as an argument."
    exit 1
fi

# Define the expected file name
fR1=$1

prefix=${fR1%_R*}

##Path and name for the results folder
out="/gpfs/gibbs/pi/braun/ERV/CM025/erv_natmed_results"

##Crosscheck whether the folder was created
if [ ! -d "$out" ]; then
	mkdir  -p $out
fi

DIR=$out/$prefix

op=$DIR/$prefix

if [ ! -d "$DIR" ]; then
	mkdir -p $DIR
fi

fR2=$prefix"_R2_001.fastq.gz"

##Path to folder with all the references
reference="/home/dp935/project/ERVtools/erv-pipe/hg19/terra"

## need to have the bai index file for bam file before sorting
## if no bai file, run the following command
##samtools index $bam
## 1. sort the bam file by query name
#samtools sort -@ ${SLURM_CPUS_PER_TASK} -n -o $op".sorted_by_name.bam" $bam

# 2. run python script, convert bam to Fastq files, 3 outputs: pair1.fq, pair2.fq & single.fq
#python /home/dp935/project/ERVtools/erv-pipe/write2fq.py $op".sorted_by_name.bam"
#echo "2. write2fq.py :"
#ls -lhtr $DIR
#echo " make sure single.fq is empty"
#echo "Make sure number of lines in single.fq is zero, otherwise input bam is not paired!"
#wc -l $DIR/single.fq
# delete intermediate large files to save some space
#rm $op".sorted_by_name.bam.*"

# 3. alignment to HUMAN TRANSCRIPTOME REFERENCE using bowtie2/2.3.4.3
# only report the best alignment using default setting (-x), and don't suppress failed alignments in sam (don't use --no-unal setting)
bowtie2 -p 10 -t -x $reference/human_index -1 $fR1 -2 $fR2 -S $DIR/alignment_human.sam
# delete intermediate large files to save some space
echo "3. aligned_to_human_transcripts :"
ls -lhtr $DIR/
# delete all human transcriptome index
rm $DIR/human_index.*.bt2

# 4. run python script, to filter out pair-end perfectly matched and pair-end 1 mismatched alignments to human transcripts reference. 2 outputs: properpair_keep.sam (keep bad alignments) & properpair_filtered.sam (filtered out pair-end pm and 1mm alignments) 
python /home/dp935/project/ERVtools/erv-pipe/filter_humanX.py $DIR/
# delete pair1.fq, pair2.fq and single.fq
echo "4. filter out humanX alignments :"
ls -lhtr $DIR/
#rm $DIR/*.fq
#rm $DIR/alignment_human.sam
# make sure no common read names between two files
samtools view $DIR/properpair_keep.sam | cut -f1 | sort | uniq > $DIR/keep_names.txt
samtools view $DIR/properpair_filtered.sam | cut -f1 | sort | uniq > $DIR/filtered_names.txt
echo "Make sure the following number of lines is zero, otherwise not filtering correctly"
comm -12 $DIR/keep_names.txt $DIR/filtered_names.txt | wc -l

# 5. extract read names of perfectly matched only on one mate from properpair_keep.sam into a text file, for final filtering out any single_pm alignments to humanX. Note: we don't filter out those single_pm alignments to humanX at this step. Because those alignments may still align to ERV-index perfectly. We just keep those reads names for future filtering.
# grep "AS:i:0" is to grep all the kept reads but still perfectly aligned to humanX
# cut -f1 is to only save read names. "sort | uniq" is to save only unique read names.
samtools view $DIR/properpair_keep.sam | grep "AS:i:0" | cut -f1 | sort | uniq > $DIR/single_pm_names.txt

# 6. sort properpair_keep.sam by name
samtools sort -n -o $DIR/alignment_merge_byname.sam $DIR/properpair_keep.sam
# delete properpair_keep.sam, properpair_human_sc.sam
echo "6. sort by name :"
ls -lhtr $DIR/
rm $DIR/properpair_human_sc.sam
# delete properpair_filtered.sam, single_filtered.sam
rm $DIR/properpair_filtered.sam
rm $DIR/properpair_keep.sam

# 7. convert the alignment_merge_byname.sam into fastq files: herv1.fq, herv2.fq, and hervS.fq (single)
python /home/dp935/project/ERVtools/erv-pipe/write2fq_human.py $DIR/
echo "7. write2fq_human.py :"
ls -lhtr $DIR/
echo "Make sure the number of lines in hervS.fq is zero!"
wc -l $DIR/hervS.fq

# 8. align herv1.fq and herv2.fq to hervquant_3173 reference, keep all MAPPED alignments
# -a: find all the alignments, --no-unal: suppress reports to sam the alignments that failed to align
bowtie2 -p 10 -t -a --no-unal -x $reference/herv_index -1 $DIR/herv1.fq -2 $DIR/herv2.fq -S $DIR/alignment.sam
echo "8. align to herv_index :"
ls -lhtr $DIR/
# delete herv1.fq, herv2.fq, hervS.fq
rm $DIR/*.fq
# delete all the index files
#rm $DIR/*_index.*.bt2
# delete all the fasta reference files
rm $DIR/*.fa


# 9. separate properly_paired alignments from only has single mate alignment
samtools view -h -f2 -o $DIR/properpair.sam -U $DIR/single.sam $DIR/alignment.sam
echo "9. separate alignment.sam :"
ls -lhtr $DIR/
rm $DIR/alignment.sam

# 10. filtering 2 conditions from pair-end and 1 condition from single-end 
python /home/dp935/project/ERVtools/erv-pipe/filter_alignments_hervquant.py $DIR/
echo "10. filter_alignments_hervquant.py :"
ls -lhtr $DIR/
rm $DIR/properpair.sam
rm $DIR/single.sam

# 11. filtering out single-end again
# combine pair_names and single_pm_names and only save unique read names
# pair_names contains the paired-end best alignments to herv_index, should be filtered out insingle_pm.sam
# single_pm_names contains the best single-end alignments to humanX, should also be filtered out in single_pm.sam
cat $DIR/pair_names.txt $DIR/single_pm_names.txt | sort | uniq > $DIR/single_filtered_names.txt
# extract the read names from single_pm_old.sam, these names need some filtering
samtools view $DIR/single_pm_old.sam | cut -f1 | sort | uniq > $DIR/single_pm_old_names.txt
# finally find the unique names in single_pm_old_names, which are not in single_filtered_names.txt
comm -23 $DIR/single_pm_old_names.txt $DIR/single_filtered_names.txt > $DIR/single_pm.txt
# sort single_pm_old.sam by name
samtools sort -n -o $DIR/single_pm_sorted.sam $DIR/single_pm_old.sam
# run python script to filter out common names between single_pm_old_names and single_filtered_names ( or only save the reads in single_pm.txt)
python /home/dp935/project/ERVtools/erv-pipe/filter_single.py $DIR/
echo "11. filter_single.py :"
ls -lhtr $DIR

############# keep duplicated version results   ################################
# 12. HTSeq counts
# pair-ends alignments perfectly matched
htseq-count --stranded=no --mode=union --nonunique=all --secondary-alignments=score --supplementary-alignments=score -a 0 -t gene -r name --samout=$DIR/pair_count.sam $DIR/properpair_pm.sam $reference/hervquant.gtf > $DIR/pair_count.txt
# pair-ends alignments only 1 mismatch
htseq-count --stranded=no --mode=union --nonunique=all --secondary-alignments=score --supplementary-alignments=score -a 0 -t gene -r name --samout=$DIR/pair1mm_count.sam $DIR/properpair_1mm.sam $reference/hervquant.gtf > $DIR/pair1mm_count.txt
# single-end alignments
htseq-count --stranded=no --mode=union --nonunique=all --secondary-alignments=score --supplementary-alignments=score -a 0 -t gene -r name --samout=$DIR/single_count.sam $DIR/single_pm.sam $reference/hervquant.gtf > $DIR/single_count.txt

echo "12. HTSeq count :"
ls -lhtr $DIR

# calculate the statistics of the flags of original bam file
samtools flagstat $DIR/alignment_merge_byname.sam > $DIR/stats.txt
rm $DIR/alignment_merge_byname.sam

# 13. count reads
python /home/dp935/project/ERVtools/erv-pipe/count_reads_hervquant.py $DIR/
#echo "13. count_reads_hervquant.py :"
ls -lhtr $DIR
rm $DIR/*_count.sam
rm $DIR/*_count.txt
# rename the 3 output files
mv $DIR/score.txt $DIR/$prefix.score.keepdup.txt
mv $DIR/score_pair.txt $DIR/$prefix.score_pair.keepdup.txt
mv $DIR/score_all.txt $DIR/$prefix.score_all.keepdup.txt

############# remove duplicated version results ################################
# 14. Change flags and remove duplicates for pair_pm, pair_1mm and single_pm
# first sort them by coordinate
samtools sort -o $DIR/properpair_pm.sortbycoord.sam $DIR/properpair_pm.sam
samtools sort -o $DIR/properpair_1mm.sortbycoord.sam $DIR/properpair_1mm.sam
samtools sort -o $DIR/single_pm.sortbycoord.sam $DIR/single_pm.sam

# then change flags: for properpair_pm and _1mm, turn off flags for "256:secondary alignment" and "2048:supplementary alignment". For single_pm, turn off all flags with "pair" or "mate" as well as secondary and supplementary flags. And then also change RNEXT="*", PNEXT="0" and TLEN="0".
python /home/dp935/project/ERVtools/erv-pipe/change_read_flags.py $DIR/
# then remove duplictes 
java -jar /home/dp935/project/ERVtools/erv-pipe/picard.jar MarkDuplicates I=$DIR/properpair_pm.sortbycoord.changeflag.sam O=$DIR/properpair_pm.rmdup.sam M=$DIR/properpair_pm.metrics.txt ASO=coordinate REMOVE_DUPLICATES=TRUE

java -jar /home/dp935/project/ERVtools/erv-pipe/picard.jar MarkDuplicates I=$DIR/properpair_1mm.sortbycoord.changeflag.sam O=$DIR/properpair_1mm.rmdup.sam M=$DIR/properpair_1mm.metrics.txt ASO=coordinate REMOVE_DUPLICATES=TRUE

java -jar /home/dp935/project/ERVtools/erv-pipe/picard.jar MarkDuplicates I=$DIR/single_pm.sortbycoord.changeflag.sam O=$DIR/single_pm.rmdup.sam M=$DIR/single_pm.metrics.txt ASO=coordinate REMOVE_DUPLICATES=TRUE
# sort the rmdup.sam by name for faster HTSeq-count and less memory buffer
# Note: I replaced the properpair_pm.sam(removed duplicates, sort by name)  with old ones (keep duplicates, sort by name)
samtools sort -n -o $DIR/properpair_pm.sam $DIR/properpair_pm.rmdup.sam
samtools sort -n -o $DIR/properpair_1mm.sam $DIR/properpair_1mm.rmdup.sam
samtools sort -n -o $DIR/single_pm.sam $DIR/single_pm.rmdup.sam

#echo "14. Remove duplicates :"
ls -lhtr $DIR
rm $DIR/*.sortbycoord.sam
rm $DIR/*.sortbycoord.changeflag.sam
rm $DIR/*.rmdup.sam

# 15. HTSeq counts, parameters I used are explained here:
# --stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature.
# --mode=union, for pair-end, if one read aligned to one gene, and its mate aligned to another gene (one case of multi-alignments), union will mark them as "ambiguous" and can be counted for both genes if --nonunique=all. However, for this case, if --mode=intersection_strict, this pair will be marked as "no_feature" and won't be counted. The htseq document is wrong.
# --nonunique=all, if "ambiguous" for multi-alignments, counted for all aligned genes.
# --secondary-alignments=score, --supplementary-alignments=score: to count for secondary/supplementary reads. The default is "ignore" instead of "score". The htseq document is wrong again! They said default is score, which is not. Actually for our case, since we already turned off all the secondary/supplementary flags, so it doesn't matter how we set these parameters.

# pair-ends alignments perfectly matched
htseq-count --stranded=no --mode=union --nonunique=all --secondary-alignments=score --supplementary-alignments=score -a 0 -t gene -r name --samout=$DIR/pair_count.sam $DIR/properpair_pm.sam $reference/hervquant.gtf > $DIR/pair_count.txt
# pair-ends alignments only 1 mismatch
htseq-count --stranded=no --mode=union --nonunique=all --secondary-alignments=score --supplementary-alignments=score -a 0 -t gene -r name --samout=$DIR/pair1mm_count.sam $DIR/properpair_1mm.sam $reference/hervquant.gtf > $DIR/pair1mm_count.txt
# single-end alignments
htseq-count --stranded=no --mode=union --nonunique=all --secondary-alignments=score --supplementary-alignments=score -a 0 -t gene -r name --samout=$DIR/single_count.sam $DIR/single_pm.sam $reference/hervquant.gtf > $DIR/single_count.txt

#echo "15. HTSeq count for remove duplicates :"
ls -lhtr $DIR

# 16. count reads
python /home/dp935/project/ERVtools/erv-pipe/count_reads_hervquant.py $DIR/
echo "16. count_reads_hervquant.py for remove duplicates:"
ls -lhtr $DIR
rm $DIR/properpair_pm.sam
rm $DIR/properpair_1mm.sam
rm $DIR/single_pm.sam

# rename the 3 output files
mv $DIR/score.txt $DIR/$prefix.score.rmdup.txt
mv $DIR/score_pair.txt $DIR/$prefix.score_pair.rmdup.txt
mv $DIR/score_all.txt $DIR/$prefix.score_all.rmdup.txt
mv $DIR/stats.txt $DIR/$prefix.stats.txt

# optional: delete intermediate files
# delete pair_count.txt, pair1mm_count.txt, single_count.txt
rm $DIR/*_count.txt
rm $DIR/*_count.sam
rm $DIR/*_names.txt
rm $DIR/single_pm_old.sam
rm $DIR/single_pm_sorted.sam
rm $DIR/*metrics.txt
rm $DIR/*.bam
echo "final lists :"
ls -lhtr $DIR

