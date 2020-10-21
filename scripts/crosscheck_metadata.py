#!/usr/bin/python

from Bio import SeqIO
import argparse
import time
import pandas as pd

today = time.strftime('%Y-%m-%d', time.gmtime())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter metadata not yet added in an existing metadata file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata1", required=True, help="FASTA file with pre-existing sequences")
    parser.add_argument("--metadata2", required=True, help="FASTA file with new sequences")
    parser.add_argument("--how", required=False, nargs=1, type=str,  default='separate', choices=['separate', 'input', 'mock'],
                        help="How the new sequences will be exported? In a 'separate' file; appended to the 'input' file, or not exported at all ('mock')?")
    args = parser.parse_args()

    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output = 'metadata_nextstrain_' + today + '.tsv'


    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_immune/nextstrain/run1_test/pre-analyses/"
    # metadata1 = path + "metadata_nextstrain.tsv"
    # metadata2 = path + "metadata_2020-09-17_08-55.tsv"
    # output = path + 'metadata_nextstrain_' + today + '.tsv'

    pd.set_option('display.max_columns', 500)

    print('\n### Scanning existing metadata...')
    # scan pre-existing metadata
    df1 = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')
    preexisting = list(set(df1['strain'].to_list()))
    df1 = df1[df1['strain'].isin(preexisting)]


    print('\n### Processing new entries...')
    # scan new metadata file
    df2 = pd.read_csv(metadata2, encoding='utf-8', sep='\t', dtype='str')
    new_entries = list(set(df2['strain'].to_list()))
    df2 = df2[df2['strain'].isin(new_entries)]
    df2 = df2[~df2['strain'].isin(preexisting)] # only new entries

    # append dataframes
    df3 = pd.concat([df1, df2], ignore_index='True')

    # export new dataframe
    df3.to_csv(output, sep='\t', index=False)


    print('\n### Exporting new file')
    if len(df1['strain'].to_list()) > 0:
        print('\n\t- A total of ' + str(len(df1['strain'].to_list())) + ' unique entries were found in the original file.')

    if len(df2['strain']) > 0:
        print('\t- A total of ' + str(len(df2['strain'].to_list())) + ' new entries were appended.\n')

    print('\t- Total entries: ' + str(len(df3['strain'].to_list())) + '\n')
