#!/usr/bin/env python
# coding: utf-8
import pysam
import pandas as pd

def filter_single_again(input_single, input_names, output_single):
    """This function iterates for each read in input_single samfile,
    if the read name is in input_names text file, 
    then write it in output_single samfile with header. """
    
    samfile = pysam.AlignmentFile(input_single, 'r')
    fpm = pysam.AlignmentFile(output_single, "wh", template=samfile)
    # load all the read names from input_names text file
    df = pd.read_csv(input_names, header=None)

    for read in samfile:
        # if read name is in df
        if df[0].str.contains(read.query_name).any() :
            #assert 'S' not in read.cigarstring, "There is soft clipping in this read!"
            fpm.write(read)
    
    samfile.close()


### need to specify the path and I/O filenames
#path = 'results/'

import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")

# input for single-end alignments
single = 'single_pm_sorted.sam'
# input single_names text file
single_names = "single_pm.txt"
# output names for filtering samfiles
output_single = 'single_pm.sam'          # single-end perfectly matched

# extract all the perfectly matched alignments from single reads
filter_single_again('{}{}'.format(path, single), '{}{}'.format(path, single_names), '{}{}'.format(path, output_single))
