#!/usr/bin/python

from Bio import SeqIO
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter newly sequenced genomes not yet added in an existing FASTA dataset of sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="FASTA file with pre-existing sequences")
    parser.add_argument("--new-genomes", required=True, help="FASTA file with new sequences")
    parser.add_argument("--max-missing", required=False, type=int,  default='30', help="Maximum percentage of Ns or gaps (int: 1-100)")
    parser.add_argument("--how", required=False, nargs=1, type=str,  default='separate', choices=['separate', 'append', 'mock'],
                        help="How the new sequences will be exported? In a 'separate' file; appended to the 'input' file, or not exported at all ('mock')?")
    args = parser.parse_args()

    dataset = args.dataset
    new_genomes = args.new_genomes
    max_gaps = args.max_missing
    how = args.how[0]
    genome_size = 29420

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_2ndWave/nextstrain/test/"
    # dataset = path + "sequences.fasta"
    # new_genomes = path + "new_sequences.fasta"
    # how = 'same'


    print('\n### Scanning existing sequences...\n')
    # scan pre-existing dataset
    preexisting = []
    duplicates = []
    too_small = []
    for entry in SeqIO.parse(open(dataset),'fasta'):
        id, seq = entry.description, str(entry.seq)
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        size = len(seq.replace('N', '').replace('-', ''))
        min_size = genome_size - int(genome_size * max_gaps/100)
        if size < min_size:
            # too_small.append(id)
            print('Partial genome (' + str(size) + ' bp): ' + id + ' is too short.')
        else:
            if id not in preexisting:
                preexisting.append(id)
            else:
                duplicates.append(id)
    # print(preexisting)
    print('Done!')


    preexist_simple = list(map(lambda strain: str(strain).replace('-', '').replace('_', ''), preexisting))

    # open output file
    outfile = ''
    if how == 'separate': # save in a separate file
        print('############################################')
        outfile = open(dataset.split('.')[0] + '_extra.fasta', 'w')
        outfile.write('')
    elif how == 'input':
        outfile = open(dataset, 'a')
        outfile.write('')
    else:
        print('\nNo output will be generated (mock run)\n')


    print('\n### Filtering and exporting new sequences...\n')
    # scan newly released genomes
    already_found = []
    new_entries = []
    outfile.write('')
    for entry in SeqIO.parse(open(new_genomes),'fasta'):
        id, seq = entry.description, str(entry.seq)
        strain = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        size = len(seq.replace('N', '').replace('-', ''))
        min_size = genome_size - int(genome_size * max_gaps/100)
        if size < min_size:
            print('Partial genome: ' + str(size) + ' bp ' + ' - ' + id + ' is too short, and was ignored...')
            too_small.append(id)
        else:
            if strain not in preexisting and strain.replace('-', '').replace('_', '') not in preexist_simple:
                print('+ ' + strain + ': new genome')
                entry = '>' + strain + '\n' + str(seq) + '\n'
                if how != 'mock':
                    outfile.write(entry)
                new_entries.append(strain)
            else:
                print('\t- ' + strain + ': this sequence was already downloaded. Skipping...')
                already_found.append(strain)

    if len(already_found) > 0:
        print('\n### The following sequences are already in the dataset:\n')
        for num, entry in enumerate(already_found):
            print('\t' + str(num + 1) + '. ' + entry)
        print('\nA total of ' + str(len(already_found)) + ' were not exported (to avoid duplicates).\n')

    if len(new_entries) > 0:
        print('\n### The following new sequecnes were exported:\n')
        for num, entry in enumerate(new_entries):
            print('\t' + str(num + 1) + '. ' + entry)
        print('\nA total of ' + str(len(new_entries)) + ' new sequences were found.\n')

    if len(duplicates) > 0:
        print('\n### WARNING: the original dataset has duplicates:\n')
        for num, entry in enumerate(duplicates):
            print('\t' + str(num + 1) + '. ' + entry)
        print('\nA total of ' + str(len(duplicates)) + ' duplicates were detected.\n')

    if len(too_small) > 0:
        print('\n### WARNING: the existing dataset contain partial genomes:\n')
        for num, entry in enumerate(too_small):
            print('\t' + str(num + 1) + '. ' + entry)
        print('\nA total of ' + str(len(too_small)) + ' partial genomes were detected and ignored.\n')
