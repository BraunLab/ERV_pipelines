#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/python
import pysam

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    Original version is from this link: https://www.biostars.org/p/306041/
    Outputs are read_dict[qname]
    If only one mate exists, read_dict[qname] = [read, None] or [None, read].
    If both mates exisit, read_dict[qname] = [read1, read2]
    """
    from collections import defaultdict
    
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair: #or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        # write the other mate into dictionary to get a pair
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            # delete the pair from the dictionary
            del read_dict[qname]
            

def soft_clip(input_samfile, softclip_samfile):
    """ This function is to check whether there is soft-clipping in the cigar string of an alignment.
    If there is soft-clipping, return True and save soft-clip reads in a 'softclip.sam' file.
    If no soft-clipping, return False.
    """
    samfile = pysam.AlignmentFile(input_samfile, 'r')
    sc = False
    for read in samfile:
        # if exists soft clip, save all soft clip reads in a samfile
        if 'S' in read.cigarstring:
            # create a new samfile
            fsc = pysam.AlignmentFile(softclip_samfile, "wh", template=samfile)
            fw.write(read)
            sc = True
    samfile.close()
    return sc


def filter_pairs(input_samfile, output_pm, output_1mm, output_names):
    """This function finds the best alignments for both ends properly mapped pair:
    If there are perfectly matched alignments for this pair, 
    save them all in output_pm samfile.
    If there isn't any perfectly matched, then look for all upto 1 mismatch alignments,
    and if there is no "N"s in the sequence, save them all in output_1mm samfile. 
    Also save all the query names (no duplicated) from output_pm and output_1mm 
    into output_names.txt file, for next step single_pm filtering."""
    
    # load the input samfile
    samfile = pysam.AlignmentFile(input_samfile, 'r')
    # create the template for writing with headers
    fpm = pysam.AlignmentFile(output_pm, "wh", template=samfile)
    f1mm = pysam.AlignmentFile(output_1mm, "wh", template=samfile)
    
    with open(output_names, mode='w') as fw:
        rp = ''     # write holder for output_pm and output_1mm
        rp2 = ''    # write holder for output_names
        for read1, read2 in read_pair_generator(samfile):
            # only look for pairs with both mates
            if (read1 is None) or (read2 is None):
                continue

            # write the pair with both ends perfectly matched into output_pm
            if read1.get_tag('AS') == 0 and read2.get_tag('AS') == 0:
                fpm.write(read1)
                fpm.write(read2)
                rp = read1.query_name
                # write the unique query name of output_pm into output_names
                if rp2 != read1.query_name:
                    fw.write(read1.query_name + "\n")
                    rp2 = read1.query_name
                    
            # if the different read pairs (means there is no pm)
            # write the pair with one end perfectly matched, the other end at most 1 mismatch
            # and no "N"s in query sequence
            elif (read1.get_tag('NM') + read2.get_tag('NM')) == 1 and \
                 (rp != read1.query_name) and \
                 ("N" not in read1.query_sequence) and ("N" not in read2.query_sequence) :
                f1mm.write(read1)
                f1mm.write(read2)
                rp = ''
                # write the unique query name of output_1mm into output_names
                if rp2 != read1.query_name:
                    fw.write(read1.query_name + "\n")
                    rp2 = read1.query_name
    
        samfile.close()
    
    
def filter_single(input_single, output_single):
    """This function finds the best alignments for single end perfectly matched alignments, 
    and write them in output_single samfile with header. """
    
    samfile = pysam.AlignmentFile(input_single, 'r')
    fpm = pysam.AlignmentFile(output_single, "wh", template=samfile)

    for read in samfile:
        # write the read with perfectly matched score
        if read.get_tag('AS') == 0 :
            #assert 'S' not in read.cigarstring, "There is soft clipping in this read!"
            fpm.write(read)
    
    samfile.close()



## need to specify the path and I/O filenames
#path = 'results/'

import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")

# input for pair-ends alignments
pair = 'properpair.sam'
# input for single-end alignments
single = 'single.sam'
# output name for softclipped checking
sc_pair = 'properpair_sc.sam'
sc_single = 'single_sc.sam'
# output names for filtering samfiles
output_pm = 'properpair_pm.sam'          # pair-end perfectly matched
output_1mm = 'properpair_1mm.sam'        # pair-end 1 mismatch
output_names = 'pair_names.txt'          # names from output_pm and output_1mm files
output_single = 'single_pm_old.sam'          # single-end perfectly matched

## check whether there is softclip in samfiles
sc1 = soft_clip('{}{}'.format(path, pair), '{}{}'.format(path, sc_pair))
sc2 = soft_clip('{}{}'.format(path, single), '{}{}'.format(path, sc_single))
print('whether there is soft-clipping reads:', sc1, sc2)

# extract all the perfectly matched (or 1 mismatched) alignments from pairs
filter_pairs('{}{}'.format(path, pair), '{}{}'.format(path, output_pm),\
             '{}{}'.format(path, output_1mm), '{}{}'.format(path, output_names))
# extract all the perfectly matched alignments from single reads
filter_single('{}{}'.format(path, single), '{}{}'.format(path, output_single))
