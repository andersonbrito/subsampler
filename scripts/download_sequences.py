# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import SeqIO
import time
import numpy as np
import argparse
import pandas as pd
import pycountry_convert as pcc
import pycountry
import os

Entrez.email = "Your.Name.Here@example.org"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Fetch newly sequenced SARS-CoV-2 genomes from NCBI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fasta", required=True, help="FASTA file with all existing genomes already downloaded")
    parser.add_argument("--metadata", required=True, help="Nextstrain metadata file")
    parser.add_argument("--skip", required=True, help="TXT file with accession number of genomes already downloaded")
    parser.add_argument("--mode", required=False, nargs=1, type=str,  default='separate', choices=['separate', 'append', 'mock'],
                        help="How the output will be exported? As separate file, or appending to existing file?")
    args = parser.parse_args()

    sequences = args.fasta
    skip = args.skip
    metadata = args.metadata
    how = args.mode[0]


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_bubble/nextstrain/run12_20201008_batch15/sampling_prop/data/'
    # sequences = path + 'gisaid_hcov-19.fasta'
    # skip = path + 'skip.txt'
    # metadata = path + 'metadata_nextstrain.tsv'
    # how = 'separate'

    # existing ncbi fasta file
    existing_ncbi = []
    for fasta in SeqIO.parse(open(sequences), 'fasta'):
        id = str(fasta.description)
        if '|' in id:
            id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        existing_ncbi.append(id)


    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    skip_list = [accno.strip() for accno in open(skip, 'r').readlines() if accno[0] not in ['\n', '#']]

    inspect = Entrez.esearch(db="nucleotide", term="txid2697049[Organism] 25000:35000[SLEN]", idtype="acc")
    total_entries = int(Entrez.read(inspect)['Count'])

    # convert state code to name
    state_names = {}
    def code2name(country, accronym):
        if country + '-' + accronym not in state_names:
            alpha2 = pcc.country_name_to_country_alpha2(country, cn_name_format="default")
            query = alpha2 + '-' + accronym
            result = pycountry.subdivisions.get(code=query).name
            state_names[country + '-' + accronym] = result
            return result
        else:
            return state_names[country + '-' + accronym]

    # open output file
    outfile = ''
    if how == 'separate': # save in a separate file
        print('############################################')
        outfile1 = open(sequences.split('.')[0] + '_ncbi.fasta', 'w')
        outfile1.write('')
        outfile2 = open(metadata.split('.')[0] + '_ncbi.tsv', 'w')
        outfile2.write('\t'.join(dfN.columns.to_list()) + '\n')

    elif how == 'append':
        outfile1 = open(sequences, 'a')
        outfile1.write('')
        outfile2 = metadata
    else:
        print('\nNo output will be generated (mock run)\n')


    ### START NCBI SEARCH

    start_at = 1
    notFound = []
    today = time.strftime('%Y-%m-%d', time.gmtime())

    # export new NCBI entries
    # output_ncbi = open(sequences.split('.')[0] + '_ncbi.fasta', 'w')
    skip_extra = open(skip, 'a')
    comment = ''

    # print(dfN['gisaid_epi_isl'][dfN['strain'] == 'England/LIVE-9B50B/2020'].values[0])

    if total_entries > 1000:
        c = 1
        for num in range(1, int(np.ceil(total_entries/1000)) + 1):
            print('\n>>> Retrieving cycle ' + str(num))
            handle = Entrez.esearch(db="nucleotide", retstart=start_at, retmax=1000,
                                    term="txid2697049[Organism] 25000:35000[SLEN]", idtype="acc")
            record = Entrez.read(handle)
            start_at = num * 1000 # define the rounds of search

            header = ['strain', 'virus', 'genbank_accession', 'date', 'country',
                      'division', 'location', 'segment', 'length', 'host', 'authors', 'date_submitted']

            previous = 1
            count = 1
            search_list = [accno.split('.')[0] for accno in record['IdList'] if accno.split('.')[0] not in skip_list]# + existing_ncbi]
            excluded = [accno.split('.')[0] for accno in record['IdList'] if accno.split('.')[0] in skip_list]# + existing_ncbi]
            c += len(excluded)

            print('A total of ' + str(len(excluded)) + ' genomes are listed in ' + skip + ', and were not re-processed in this cycle.')
            for accno in search_list:#[245:250]:
                print('\n' + str(c) + '/' + str(total_entries))

                # exporting accession numbers of processed GISAID-NCBI duplicates
                if accno in skip_list:# + existing_ncbi:  # and strain in dfN['strain'].to_list():
                    print("\t- " + accno + ': Genome already present in dataset. Skipping...')
                else:
                    accno = accno.split('.')[0]
                    new_row = {} # to update metadata dataframe

                    for column_name in header:
                        new_row[column_name] = ''

                    current = time.time()
                    delay = 0.4
                    elapsed = current - previous
                    wait = delay - elapsed

                    try:
                        # print(accno)
                        handle = Entrez.efetch(db="nucleotide", id=accno, rettype="gb", retmode="text")
                        strain, virus, genbank_accession, date, country, division,\
                        location, segment, length, host, authors, date_submitted = ['' for x in header]
                        sequence = ''
                        isolate = ''

                        for seq_record in SeqIO.parse(handle, "gb"):
                            sequence = str(seq_record.seq) # genome sequence

                            for feature in seq_record.annotations['references']:
                                length = str(len(seq_record.seq))
                                authors = feature.authors.split(",")[0] + " et al"
                                if feature.title =='Direct Submission':
                                    date_submitted = pd.to_datetime(feature.journal.split('(')[1].split(')')[0]).strftime('%Y-%m-%d')
                                    # print(date_submitted)
                            for feature in seq_record.features:
                                if feature.type == 'source':
                                    try:
                                        date = feature.qualifiers['collection_date'][0]
                                        # print(date)
                                        if len(date) > 7:
                                            date = pd.to_datetime(date).strftime('%Y-%m-%d')
                                        elif len(date) == 4:
                                            date = date + '-XX-XX'
                                        elif len(date) == 7:
                                            date = date + '-XX'
                                    except:
                                        pass
                                    try:
                                        origin = feature.qualifiers['country'][0]
                                        if ':' in origin:
                                            country = origin.split(":")[0]
                                            if len(origin.split(":")) > 0:  # get subnational location information
                                                subnational = origin.split(":")[1]
                                                if ',' in subnational:
                                                    division = subnational.split(',')[0].strip()
                                                    location = subnational.split(',')[1].strip()
                                                else:
                                                    division = subnational.strip()
                                        else:
                                            country = origin
                                    except:
                                        pass
                                    try:
                                        if len(division) == 2:
                                            division = code2name(country, division)
                                    except:
                                        pass
                                    # get strain name, unique ID
                                    try:
                                        isolate = feature.qualifiers['isolate'][0]
                                        if '/' in isolate:
                                            strain = '/'.join([country] + isolate.split('/')[-2:]).replace(' ','')
                                        else:
                                            strain = '/'.join([country] + [isolate] + [str(2020)]).replace(' ','')
                                    except:
                                        break
                                    try:
                                        host = feature.qualifiers['host'][0]
                                    except:
                                        pass

                        data = {'strain': strain, 'virus': 'ncov', 'genbank_accession': accno, 'date': date,
                                        'country': country, 'division': division, 'location': location, 'segment': segment,
                                        'length': length, 'host': host, 'authors': authors, 'date_submitted': date_submitted}
                        # print(data)
                        new_row.update(data)

                        # report missing data
                        if '' in [data['date'], data['country'], isolate]:
                            # print([data['date'], data['country'], isolate])
                            # print(accno + ': no strain, date or country information')
                            c += 1
                            if comment == '':
                                comment = '\n# Processed on ' + today + '\n'
                                skip_extra.write(comment)
                            skip_extra.write(accno + '\n')
                            skip_extra. flush()
                            print('\t- ' + accno + ': Missing strain, country or date metadata. Skipping... ')

                            continue

                        ncbi_header = 'hCoV-19/' + strain + '|' + accno + '|' + date
                        try:
                            epi_isl = dfN['gisaid_epi_isl'][dfN['strain'] == strain].values[0]
                        except:
                            epi_isl = ''


                        # exporting new metadata lines
                        if strain not in dfN['strain'].to_list():# and strain.replace('_', '-') not in dfN['strain'].to_list():
                            # exporting new NCBI sequences
                            if accno not in existing_ncbi and epi_isl in ['?', '', np.nan]:
                                if how == 'separate':
                                    outfile1.write('>' + ncbi_header + '\n' + sequence + '\n')
                                    outfile1.flush()

                                elif how == 'append':
                                    outfile1.write('>' + ncbi_header + '\n' + sequence + '\n')
                                    outfile1.flush()
                                existing_ncbi.append(accno)
                                print('\t- Exporting NCBI fasta: ' + ncbi_header)


                            dfG = pd.DataFrame(columns=dfN.columns.to_list())
                            dfG = dfG.append(new_row, ignore_index=True)
                            if how == 'separate':
                                # if os.path.getsize(outfile2) > 0:
                                #     dfG.to_csv(outfile2, sep='\t', index=False, header=False, mode='a')
                                # else:
                                dfG.to_csv(outfile2, sep='\t', index=False, header=False, mode='w')
                            elif how == 'append':
                                dfG.to_csv(outfile2, sep='\t', index=False, header=False, mode='a')

                            # print(list(dfN.loc[dfN['strain'] == strain].values))
                            print("\t- Exporting NCBI metadata: " + ncbi_header)


                        else:
                            # if not epi_isl in ['?', '', np.nan] and accno not in skip_list:
                            if comment == '':
                                comment = '\n# Processed on ' + today + '\n'
                                skip_extra.write(comment)
                            skip_extra.write(accno + '\n')
                            skip_extra.flush()
                            print('\t- ' + accno+ ': GISAID entry. Skipping... ')
                    except:
                        print("\t- " + accno + ": entry not found on NCBI, or backend searching failed")
                        notFound.append(accno)
                    previous = time.time()
                    count += 1
                    c += 1

    # list entries not found on NCBI
    if len(notFound) > 0:
        print('\nThe following genomes were not retrieved:\n')
        for entry in notFound:
            print('\t' + entry)



