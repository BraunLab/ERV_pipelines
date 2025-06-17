#!/usr/bin/env python
# coding: utf-8

import pysam

def decimalToBinary(n):
    """change an integer from decimal to a binary string."""
    return bin(n).replace("0b","")

def change_flags(flag, index=[0,3,4,5,6,8,10,11]):
    """
    Turn off flag (change from 1 to 0) at each index.
    The default is to turn off all flags related to 
    'pair', 'mate', 'secondary' and 'supplementary'.
    Index for each flag is in reverse order, as shown below:
    index 0: supplementary alignment;
    index 3: secondary alignment;
    index 4: second in pair;
    index 5: first in pair;
    index 6: mate reverse strand;
    index 8: mate unmapped;
    index 10: read mapped in proper pair;
    index 11: read paired.
    """
    
    assert type(flag) is int, "flag should be an integer!"
    assert type(index) is list, "index should be a list of integers!"
    assert set(index).issubset(set([0,1,2,3,4,5,6,7,8,9,10,11])), \
     "Your flag index is out of range, should be integers from 0 to 11!"
    assert len(index) > 0, "index list shouldn't be empty!"
    # convert flag from decimal integer to Binary string
    bs = decimalToBinary(flag)
    # pad with zeros on the left to 12 digits
    bs_long = bs.zfill(12)
    # convert padded string to list
    bsl = list(bs_long)
    # turn off 1 at position 0, 3, 4, 5, 6, 8, 10, 11 to zero
    for i in index:
        bsl[i]='0'
    # join list to string
    bsj = ''.join(bsl)
    # convert from binary string to decimal
    flag_new = int(bsj, base=2)
    return flag_new



def change_read_flags(bam_in, bam_out, index=[0,3,4,5,6,8,10,11]):
    """Change the flag of each read in bam_in, turn off the flag at each index from 1 to 0. 
    And save the new reads in bam_out. 
    Index for each flag is in reverse order, as shown below:
    index 0: supplementary alignment;
    index 3: secondary alignment;
    index 4: second in pair;
    index 5: first in pair;
    index 6: mate reverse strand;
    index 8: mate unmapped;
    index 10: read mapped in proper pair;
    index 11: read paired.
    """
    
    # read in samfile, if for bamfile, use "rb" instead of "r"
    bamfile = pysam.AlignmentFile(bam_in, 'r')
    # write into samfile template, if write into bamfile, use "wb" 
    # instead of "wh".
    fw = pysam.AlignmentFile(bam_out, "wh", template=bamfile)
    for read in bamfile:
        newflag = change_flags(read.flag, index)
        read.flag = newflag
        fw.write(read)
    bamfile.close()   



def change_single_read(bam_in, bam_out, index=[0,3,4,5,6,8,10,11]):
    """Change the flag of each read in bam_in, 
    turn off the flag at each index from 1 to 0. 
    And save the new reads in bam_out. 
    Index for each flag is in reverse order, as shown below:
    index 0: supplementary alignment;
    index 3: secondary alignment;
    index 4: second in pair;
    index 5: first in pair;
    index 6: mate reverse strand;
    index 8: mate unmapped;
    index 10: read mapped in proper pair;
    index 11: read paired.
    
    Also change column 7(RNEXT), 8(PNEXT) and 9(TLEN) to "*", "0", and "0"
    """
    # read in samfile, if read in bam, use "rb" instead of "r".
    bamfile = pysam.AlignmentFile(bam_in, 'r')
    # write into samfile template, if write into bamfile, use "wb"
    # instead of "wh".
    fw = pysam.AlignmentFile(bam_out, "wh", template=bamfile)
    for read in bamfile:
        newflag = change_flags(read.flag, index)
        # change flag
        read.flag = newflag
        # change RNEXT to "*", pysam use -1 to indicate "*"
        read.next_reference_id = -1
        # change PNEXT to "0", pysam use -1 to indicate "0"
        read.next_reference_start = -1
        # change TLEN to "0"
        read.template_length = 0
        # write new read in bam_out
        fw.write(read)
    bamfile.close()   



# change read for single_pm.sam
#path="results/"

import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")


single_in="single_pm.sortbycoord.sam"
single_out="single_pm.sortbycoord.changeflag.sam"
change_single_read("{}{}".format(path, single_in), "{}{}".format(path, single_out))

# for properpair_pm.sam and _1mm.sam, only turn off secondary and supplementary flags
pairpm_in="properpair_pm.sortbycoord.sam"
pairpm_out="properpair_pm.sortbycoord.changeflag.sam"
change_read_flags("{}{}".format(path, pairpm_in), "{}{}".format(path, pairpm_out), [0, 3])

pair1mm_in="properpair_1mm.sortbycoord.sam"
pair1mm_out="properpair_1mm.sortbycoord.changeflag.sam"
change_read_flags("{}{}".format(path, pair1mm_in), "{}{}".format(path, pair1mm_out), [0, 3])
