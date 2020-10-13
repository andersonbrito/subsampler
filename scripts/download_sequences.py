# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import SeqIO
import time
import numpy as np
import argparse
import pandas as pd
import pycountry_convert as pcc
import pycountry

Entrez.email = "Your.Name.Here@example.org"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Append newly sequenced genomes to current genome dataset, and export metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fasta", required=True, help="FASTA file with all existing genomes already downloaded")
    parser.add_argument("--skip", required=True, help="TXT file with accession number of genomes already downloaded")
    parser.add_argument("--metadata", required=True, help="Newly generated metadata file")
    args = parser.parse_args()

    ncbi_fasta = args.fasta
    redundant = args.skip
    metadata = args.metadata

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_2ndWave/nextstrain/test/'
    # ncbi_fasta = path + 'sequences.fasta'
    # redundant = path + 'skip.txt'
    # metadata = path + 'metadata_short.tsv'


    # existing ncbi fasta file
    existing_ncbi = []
    for fasta in SeqIO.parse(open(ncbi_fasta), 'fasta'):
        id = str(fasta.description).split('|')[1]
        existing_ncbi.append(id)

    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    # ns_entries = dfN['strain'].to_list()
    dup_seqs = [accno.strip() for accno in open(redundant, 'r').readlines() if accno[0] not in ['\n', '#']]

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

    ### START NCBI SEARCH

    start_at = 1
    notFound = []
    today = time.strftime('%Y-%m-%d', time.gmtime())

    # export new NCBI entries
    export_ncbi = open(ncbi_fasta, 'a')
    new_redundants = open(redundant, 'a')
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
            search_list = [accno.split('.')[0] for accno in record['IdList'] if accno.split('.')[0] not in dup_seqs + existing_ncbi]
            excluded = [accno.split('.')[0] for accno in record['IdList'] if accno.split('.')[0] in dup_seqs + existing_ncbi]
            c += len(excluded)

            print('A total of ' + str(len(excluded)) + ' genomes are listed in ' + redundant + ', and were not re-processed in this cycle.')
            for accno in search_list:#[245:250]:
                # print(accno)
                print('\n' + str(c) + '/' + str(total_entries))

                # exporting accession numbers of processed GISAID-NCBI duplicates
                if accno in dup_seqs + existing_ncbi:  # and strain in dfN['strain'].to_list():
                    print("\t- " + accno +  ': Genome already processed')
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
                                    print(date_submitted)
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
                                new_redundants.write(comment)
                            new_redundants.write(accno + '\n')
                            new_redundants. flush()
                            print('\t- ' + accno + ': Missing strain, country or date metadata. Skipping... ')

                            continue

                        ncbi_header = 'hCoV-19/' + strain + '|' + accno + '|' + date
                        try:
                            epi_isl = dfN['gisaid_epi_isl'][dfN['strain'] == strain].values[0]
                        except:
                            epi_isl = ''

                        # # exporting accession numbers of processed GISAID-NCBI duplicates
                        # if accno in dup_seqs + existing_ncbi:# and strain in dfN['strain'].to_list():
                        #         print(str(c) + '/' + str(total_entries) + " - Sequence already processed: " + accno)

                        # exporting new metadata lines
                        if strain not in dfN['strain'].to_list():# and strain.replace('_', '-') not in dfN['strain'].to_list():
                            dfG = pd.DataFrame(columns=dfN.columns.to_list())
                            dfG = dfG.append(new_row, ignore_index=True)
                            dfG.to_csv(metadata, sep='\t', index=False, header=False, mode='a')

                            # print(list(dfN.loc[dfN['strain'] == strain].values))
                            print("\t- Exporting NCBI metadata: " + ncbi_header)

                            # exporting new NCBI sequences
                            if accno not in existing_ncbi and epi_isl in ['?', '', np.nan]:
                                export_ncbi.write('>' + ncbi_header + '\n' + sequence + '\n')
                                export_ncbi.flush()
                                existing_ncbi.append(accno)
                                print('\t- Exporting NCBI fasta: ' + ncbi_header)

                        else:
                            # if not epi_isl in ['?', '', np.nan] and accno not in dup_seqs:
                            if comment == '':
                                comment = '\n# Processed on ' + today + '\n'
                                new_redundants.write(comment)
                            new_redundants.write(accno + '\n')
                            new_redundants.flush()
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



