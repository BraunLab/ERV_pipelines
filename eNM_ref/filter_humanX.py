#!/usr/bin/env python
# coding: utf-8

import pysam

def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found, 
    including both proper_paired and unproper_paired. 
    Original version is from this link: https://www.biostars.org/p/306041/
    Outputs are read_dict[qname]
    If only one mate exists, read_dict[qname] = [read, None] or [None, read].
    If both mates exisit, read_dict[qname] = [read1, read2]
    """
    from collections import defaultdict
    
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        # select all pairs, including both proper_paired and unproper_paired
        if not read.is_paired: 
        # if want to select only proper_paired, use the following if command
        #if not read.is_proper_pair: # or read.is_secondary or read.is_supplementary:
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
    # create a new samfile fsc for writing soft-clipped sequences
    fsc = pysam.AlignmentFile(softclip_samfile, "wh", template=samfile)
    
    sc = False
    for read in samfile:
        # if exists soft clip, save all soft clip reads in a samfile
        if (read.cigarstring != None) and ('S' in read.cigarstring):
            fsc.write(read)
            sc = True
    samfile.close()
    return sc




def filterout_pairs(input_samfile, filtered_pair, keep_pair):
    """This function excludes the best alignments for both ends properly mapped pair:
    If there are perfectly matched alignments for this pair, 
    remove them all from keep_pair samfile, and write them into filtered_pair samfile.
    If there isn't any perfectly matched, then look for all upto 1 mismatch alignments, 
    remove them all from keep_pair samfile, and write them into filtered_pair samfile. 
    Then for the bad alignments, write them into keep_pair samfile."""
    
    # load the input samfile
    samfile = pysam.AlignmentFile(input_samfile, 'r')
    # create the template samfile for writing bad alignments with headers
    fk = pysam.AlignmentFile(keep_pair, "wh", template=samfile)
    # create the template samfile for writing filtered out pairs with headers
    ff = pysam.AlignmentFile(filtered_pair, "wh", template=samfile)
    
    for read1, read2 in read_pair_generator(samfile):
        # only look for pairs with both mates
        if (read1 is None) or (read2 is None):
            continue

        # exclude the pair with both ends perfectly matched or with 1 mismatch
        # 1 mismatch means 1 edit distance to reference
        if (read1.has_tag('AS') and read2.has_tag('AS') and \
            read1.get_tag('AS') == 0 and read2.get_tag('AS') == 0) or \
           (read1.has_tag('NM') and read2.has_tag('NM') and \
           (read1.get_tag('NM') + read2.get_tag('NM')) == 1):
            # write the filtered out (best alignments) pairs into filtered_pair samfile
            ff.write(read1)
            ff.write(read2)
        else:
            # write the unfiltered pairs (bad alignments) into keep_pair samfile
            fk.write(read1)
            fk.write(read2)

    samfile.close()


# didn't use this filterout_single function
def filterout_single(input_single, output_single):
    """This function excludes the best alignments for single end perfectly matched alignments, 
    and remove them from output_single samfile with header. """
    
    samfile = pysam.AlignmentFile(input_single, 'r')
    fpm = pysam.AlignmentFile(output_single, "wh", template=samfile)

    for read in samfile:
        # write the read DOES NOT have perfectly matched score
        if (not read.has_tag('AS')):
            fpm.write(read)
        elif read.get_tag('AS') != 0 :
            fpm.write(read)
    
    samfile.close()


# didn't use this create_gtf function
def create_gtf(input_fasta, output_gtf, source = 'GenBank', feature = 'gene', start = '1',                score = '.', strand = '+', frame = '0'):
    """This function create a GTF file from the indexed FASTA file."""

    ### default settings for GTF features or rows:
    #source = 'GenBank'
    #feature = 'gene'
    #start = '1'  # Start position of the feature, with sequence numbering starting at 1.
    #score = '.'  # a floating point value
    #strand = '+' # defined as + (forward) or - (reverse).
    ## One of '0', '1' or '2'. '0' indicates that the first base of the feature 
    ## is the first base of a codon, '1' that the second base is the first base of a codon, 
    ## and so on..
    #frame = '0'   
    ## A semicolon-separated list of tag-value pairs, 
    ## providing additional information about each feature.
    ##attribute = '.'
             
    with open(input_fasta, 'r') as f, open(output_gtf, mode='w') as fw:
        n = 1
        for line in f:
            # the header line for each gene or erv, eg: '>X57147|ERV9-1\n'
            if line[0] == '>':
                # gene name, remove the end '\n'
                seqname = line[1:-1]
            # the genome sequence line for each gene or erv
            else:
                # End position for the gene, with sequence numbering starting at 1.
                # line ends with '\n', so .rstrip() to remove the '\n'
                end = str(len(line.rstrip()))
                # attribute is like: gene_id "100";
                attribute = 'gene_id' + ' ' + '"' + str(n) + '"' + ';'
                # gtf content contains 9 columns: 
                # each column is seperated by '\t', and each feature(or row) is end with '\n'
                content = [seqname, source, feature, start, end, \
                           score, strand, frame, attribute, '\n']
                fw.write('\t'.join(content))
                # update n for gene_id
                n += 1
    # after finished writing, print the last row to check
    print('Last row in GTF:')
    print('\t'.join(content))




################ need to specify the path and I/O filenames
import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")

#path = 'results/'
# input for samfile for filtering
pair = 'alignment_human.sam'
#pair = 'properpair_human.sam'  # use this if your single.fq is not empty
## input for single-end alignments
#single = 'single_human.sam'
# output name for softclipped checking
sc_pair = 'properpair_human_sc.sam'
#sc_single = 'single_human_sc.sam'
# output names for filtered samfiles
filtered_pair = 'properpair_filtered.sam'    # pair-end perfectly matched and 1 mismatched 
# output names for keep bad alignments samfiles
keep_pair = 'properpair_keep.sam'        # pair-end perfectly matched and 1 mismatched out

## check whether there is softclip in samfiles
sc1 = soft_clip('{}{}'.format(path, pair), '{}{}'.format(path, sc_pair))
#sc2 = soft_clip('{}{}'.format(path, single), '{}{}'.format(path, sc_single))
print('whether there is soft-clipping reads:', sc1)

# exclude all the perfectly matched and 1 mismatched alignments from pairs
filterout_pairs('{}{}'.format(path, pair), '{}{}'.format(path, filtered_pair), '{}{}'.format(path, keep_pair))
## exclude all the perfectly matched alignments from single reads
#filterout_single('{}{}'.format(path, single), '{}{}'.format(path, output_single))
