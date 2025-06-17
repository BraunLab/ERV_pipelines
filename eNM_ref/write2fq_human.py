#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/python

import pysam

def reverseComplement(s):
    '''return the reverse complement string rc of a sequence s.'''
    assert type(s) is str, "Input variable should be a string!"
    # define the complement dictionary
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    # initialize the reverse complement string as empty string
    rc = ''
    # if the input string has lower cases, convert them to upper case
    for base in s.upper():
        assert base in complement, "Your sequence contains other characters except ATCGN"
        # first complement, then reverse
        rc = complement[base] + rc   # the order of 'complement + rc' is very important
    return rc



def QtoPhred33(Q):
    """Given a list of quality scores as integer, return a string of quality in Phred33 format.
    Note: This Q is in the format of L - Illumina 1.8+ Phred+33,  raw reads typically [0, 41]."""
    # assert type(Q) is list, "Your input variable should be a list!"
    assert type(Q[0]) is int, "Each item in your input list should be an integer!"
    assert min(Q) >= 0 and max(Q) <= 41, "Your quality intergers should be in the range [0, 41]!"
    return ''.join([chr(q+33) for q in Q])



def reverse(string): 
    """Return the reverse order of a string. Doesn't change the original string. """
    string = string[::-1] 
    return string



def compareID(read1, read2):
    """
    This function needs the pysam module. Make sure you installed it following this link:
    https://pysam.readthedocs.io/en/latest/installation.html
    
    This function is to compare the IDs between read1 and read2 and return 3 different status:
        Return False: If two reads have different ID, means they are not paired;
        Return True: If two reads have same ID but different tags, means they are paired;
        Return None: If two reads have same ID and same tags, means they are same read.
    read1 and read2 are in pysam.AlignedSegment format.
    read1.query_name is like 'UNC12-SN629_87:4:1101:1150:10255/1' or "UNC12-SN629_87:4:1101:1150:10255"
    
    This function first check whether the ID has tags "/1" and "/2".
    If has tags, use tags to compare; If doesn't have tags, use is.read1 and is.read2 to compare.
    """
    
    if read1.query_name[-2]=='/': # IDs have tags
        # different IDs, return False
        if read1.query_name[:-1] != read2.query_name[:-1] : return False
        # same ID different tags, same pair different mates, return True
        elif read1.query_name[-1] != read2.query_name[-1] : return True
        # same ID and same tags, same reads, return None
        else: return None
    else: # IDs don't have tags
        # different IDs, return False
        if read1.query_name != read2.query_name : return False
        # same ID different tags, same pair different mates, return True
        elif (read1.is_read1 and read2.is_read2) or (read1.is_read2 and read2.is_read1): return True
        # same ID and same tags, same reads, return None
        else: return None



def write_fastq(read, fw):
    """
    This function needs the pysam module. Make sure you installed it following this link:
    https://pysam.readthedocs.io/en/latest/installation.html
    
    This function is to write a 'read' (in pysam.AlignedSegment format) into a fastq file.
    'fw' is the writing fastq file handle. 'fw' must be defined before running this function.
    An example of how to define 'fw' is like the following command:
    with open('write_fastq_file_name.fq', mode='w') as fw:
    """
    
    ## Check if read is on reverse strand
    # The pysam.AlignedSegment class has many useful features, 
    # please refer 'pysam documentation, Release 0.15.0, page 12'.   
    # '.is_reverse': True if read is mapped to reverse strand.
    if read.is_reverse:
        # write read name. '.query_name': the query template name (None if not present).
        fw.write('@' + read.query_name + "\n")                     
        # write reverse complement of sequence
        fw.write(reverseComplement(read.query_sequence) + "\n")   
        fw.write('+' + "\n")
        # write reverse of quality, quality is transferred to Phred33 format
        fw.write(reverse(QtoPhred33(read.query_qualities)) + "\n") 
                
    # if read is not on reverse strand, write the same sequence and quality
    else: 
        fw.write('@' + read.query_name + "\n")
        fw.write(read.query_sequence + "\n")                      
        fw.write('+' + "\n")
        fw.write(QtoPhred33(read.query_qualities) + "\n")



def write2fastq(sam_filename, pair1_fastq_filename, pair2_fastq_filename, single_fastq_filename):
    assert type(sam_filename) is str
    assert type(pair1_fastq_filename) is str
    assert type(pair2_fastq_filename) is str
    assert type(single_fastq_filename) is str
    # read the sam file ('rb', -r: read, -b:bam, -u:uncompressed, -w:write, -h:output header info)
    # Note: make sure the bam/sam file is sorted by name!
    samfile = pysam.AlignmentFile(sam_filename, 'r')
    # write into 3 fastq files
    with open(pair1_fastq_filename, mode='w') as fp1,          open(pair2_fastq_filename, mode='w') as fp2,          open(single_fastq_filename, mode='w') as fs:
        # initialize previous reads holder 'rp'
        rp1 = []
        rp2 = []
        wp = False
        ws = False
        for read in samfile:
            if not rp1: rp1 = read
            elif not rp2: 
                rp2 = read
                if compareID(rp1,rp2) == None: rp2 = []
                elif compareID(rp1,rp2): wp = True  
                else: ws = True
            else:
                if compareID(rp2, read) == None: continue
                elif compareID(rp2, read) and ws:
                    wp = True
                    write_fastq(rp1, fs)
                    ws = False
                    rp1 = rp2
                    rp2 = read
            # be careful with this, might have error
                else:
                    if wp:
                        if rp1.is_read1 and rp2.is_read2:
                            write_fastq(rp1, fp1)
                            write_fastq(rp2, fp2) 
                        if rp1.is_read2 and rp2.is_read1:
                            write_fastq(rp2, fp1)
                            write_fastq(rp1, fp2)
                        wp = False
                        rp1 = read
                        rp2 = []
                    if ws:
                        write_fastq(rp1, fs)
                        write_fastq(rp2, fs)
                        ws = False
                        rp1 = read
                        rp2 = []
        
        if rp2:
            if compareID(rp1,rp2) == None: write_fastq(rp1, fs)
            elif compareID(rp1, rp2):
                if rp1.is_read1 and rp2.is_read2:
                    write_fastq(rp1, fp1)
                    write_fastq(rp2, fp2)
                if rp1.is_read2 and rp2.is_read1:
                    write_fastq(rp2, fp1)
                    write_fastq(rp1, fp2)
            else:
                write_fastq(rp1, fs)
                write_fastq(rp2, fs)
        else: write_fastq(rp1, fs)
            
    samfile.close()



# load the input bam file from ./results folder
# then write the two mates of each pair into pair1.fq and pair2.fq 
# if there exists single read, (without mate), write it into single.fq

# e.g. filename = './results/UNCID_1187919'
#path = 'results/'

import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")


write2fastq('{}alignment_merge_byname.sam'.format(path), '{}herv1.fq'.format(path), '{}herv2.fq'.format(path), '{}hervS.fq'.format(path))
