#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/python

def extract_total_reads(stats_file):
    """
    This function extracts the first number from the stats_file, and return with that number.
    """
    import re
    with open(stats_file, 'r') as f:
        n = 0
        # only read the first line
        for line in f:
            if n > 0:
                break
            # extract the first number in line
            totals = re.match(r'(\d+)', line)
            # convert it to float type for further normalization
            totals = float(totals.group())
            n += 1
    return totals

######################## need to specify I/O filenames
#path = 'results/'

import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")


stats_file = '{}stats.txt'.format(path)
tr = extract_total_reads(stats_file)
print('total number of reads is:', tr)



def count_reads(total_reads, output_score, gtf_file, pair_pm, pair_1mm = None, single = None, p = 8):
    """
    This function is to count the number of aligned reads on each gene, 
    and generate the 'output_score' file. 
    There are 3 possible input files for aligned reads:
    'pair_pm' is required, perfectly matched pair-ends reads.
    'pair_1mm' is optional, only one mismatched pair-ends reads.
    'single' is optional, perfectly matched single-end reads.
    The input 'gtf_file' is required, which is the reference gene (erv) file.
    
    The score is calculated as:
    score = (counts of pair_pm * 2) + (counts of pair_1mm * 2) + (counts of single)
    normal = (score/total_reads) * (10 ** p)
    
    The output_score file has 3 columns: 'gene_name', 'gene_id' and 'counts'.
    3 new output text files will be created based on 3 input read_counts files.
    e.g. A new pair_pm file will called 'pair_pm_new.txt', which has the gene_name column added.
    
    p: power to show the number in 1-100 range, default is 8. 
    """
    import pandas as pd
    
    if type(total_reads) is not float: total_reads=float(total_reads)
    
    ## load GTF file, only load column 0 and 8, and name them gene_name, gene_id
    df1 = pd.read_table(gtf_file, usecols=[0,8], names=('gene_name','gene_id'))
    # extract the number from 'gene_id' column
    df1['gene_id'] = df1['gene_id'].str.extract(r'(\d+)')
    # number of total genes (rows) in GTF file
    n = df1.shape[0]
        
    ## load pair_pm file, without header, rename the columns
    df2 = pd.read_table(pair_pm, names=('gene_id', 'counts'))
    # outer join the two on gene_id
    df2new = df1.merge(df2, how='outer')
    # write it into pair_pm_new file
#    df2new.to_csv(pair_pm[:-4] + '_new.txt', sep='\t', index=False)
    # count the reads
    df_score = df2new[:n]
#    df_score[['counts']] = df_score[['counts']]*2
    # Fixed code:
    df_score.loc[:, 'counts'] = df_score['counts']*2
#    df_score[['normal']] = (df_score[['counts']]/total_reads) * (10**p)
    
    # Fixed code:
    df_score.loc[:, 'normal'] = (df_score['counts']/total_reads) * (10**p)

    # if has pair_1mm:
    if pair_1mm:
        df3 = pd.read_table(pair_1mm, names=('gene_id', 'counts'))
        # outer join the two on gene_id
        df3new = df1.merge(df3, how='outer')
        # write it into pair_1mm file
#        df3new.to_csv(pair_1mm[:-4] + '_new.txt', sep='\t', index=False)
        # count the reads
        df3_score = df3new[:n]
#        df_score[['counts']] = df_score[['counts']] + df3_score[['counts']] * 2
    # Fixed code:
        df_score.loc[:, 'counts'] = df_score['counts'] + df3_score['counts'] * 2        
        df_score[['normal']] = (df_score[['counts']]/total_reads) * (10**p)
        
    # if has single:
    if single:
        df4 = pd.read_table(single, names=('gene_id', 'counts'))
        # outer join the two on gene_id
        df4new = df1.merge(df4, how='outer')
        # write it into pair_pm file
#        df4new.to_csv(single[:-4] + '_new.txt', sep='\t', index=False)
        # count the reads
        df4_score = df4new[:n]
        df_score[['counts']] = df_score[['counts']] + df4_score[['counts']]
        df_score[['normal']] = (df_score[['counts']]/total_reads) * (10**p)
    
    # finally write the df_score into text file
    df_score.to_csv(output_score, sep='\t', index=False)



###################### need to specify I/O filenames
#path = 'results/'

import sys

if len(sys.argv) > 1:
    path = sys.argv[1]  # Take the first argument as the path
    print(f"The provided path is: {path}")
else:
    print("No path provided. Please run the script with a path argument.")

output_score1 = '{}score_all.txt'.format(path)
output_score2 = '{}score.txt'.format(path)
output_score3 = '{}score_pair.txt'.format(path)
gtf_file = '/home/dp935/project/ERVtools/erv-pipe/hg19/hervquant.gtf'
pair_pm = '{}pair_count.txt'.format(path)
pair_1mm = '{}pair1mm_count.txt'.format(path)
single = '{}single_count.txt'.format(path)

count_reads(tr, output_score1, gtf_file, pair_pm, pair_1mm, single)
count_reads(tr, output_score2, gtf_file, pair_pm, pair_1mm = None, single = None)
count_reads(tr, output_score3, gtf_file, pair_pm, pair_1mm, single = None)
