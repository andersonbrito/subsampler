import pandas as pd
from Bio import SeqIO
from epiweeks import Week
import random
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file with genomes named as in metadata")
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--genome-matrix", required=True, help="TSV file showing corrected genome counts per epiweek")
    parser.add_argument("--max-missing", required=True, type=int,  help="Maximum percentage of Ns or gaps (int = 1-100)")
    parser.add_argument("--refgenome-size", required=True, type=int,  help="Reference genome size")
    parser.add_argument("--keep", required=False, help="List of samples to keep, in all instances")
    parser.add_argument("--remove", required=False, help="List of samples to remove, in all instances")
    parser.add_argument("--seed", required=False, type=int,  help="Seed number for pseudorandom sampling of genomes")
    parser.add_argument("--index-column", required=True, help="Metadata column with unique geographic information")
    parser.add_argument("--date-column", required=True, type=str,  help="Metadata column containing the collection dates")
    parser.add_argument("--filter-column", required=False, type=str,  help="Column with dates used to filter by date")
    parser.add_argument("--start-date", required=False, type=str,  help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str,  help="End date in YYYY-MM-DD format")
    parser.add_argument("--sampled-sequences", required=True, help="Sampled genomes")
    parser.add_argument("--sampled-metadata", required=True, help="Sampled metadata")
    parser.add_argument("--report", required=True, help="List of statistics related to the sampling scheme")
    args = parser.parse_args()


    input1 = args.sequences
    input2 = args.metadata
    input3 = args.genome_matrix
    genome_size = args.refgenome_size
    max_gaps = args.max_missing
    keep = args.keep
    remove = args.remove
    seed = args.seed
    geo_level = args.index_column
    date_col = args.date_column
    filter_col = args.filter_column
    start_date = args.start_date
    end_date = args.end_date
    output1 = args.sampled_sequences
    output2 = args.sampled_metadata
    output3 = args.report


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/github/subsampler/'
    # input1 = path + 'data/sequences_corrected.fasta'
    # input2 = path + 'data/metadata.tsv'
    # input3 = path + 'outputs/matrix_cases_epiweeks_corrected.tsv'
    # keep = path + 'config/keep.txt'
    # remove = path + 'config/remove.txt'
    # geo_level = 'iso'
    # seed = 2007
    # # seed = None
    # date_col = 'date'
    # filter_col = 'date'
    # genome_size = 29420
    # max_gaps = 5
    # start_date = '2019-12-15'
    # end_date = '2020-08-29'
    # # start_date = None
    # # end_date = None
    # output1 = path + 'sequences_corrected0001.fasta'
    # output2 = path + 'metadata_corrected0001.tsv'
    # output3 = path + 'report.txt'


    if seed == None:
        seed = random.random()

    separator = '\t'
    if seed == None:
        seed = random.random()

    if filter_col == None:
        filter_col = date_col

    # pd.set_option('display.max_columns', 500)

    print('\n### Loading sequences...')
    # get sequence headers
    fasta_headers = []
    for fasta in SeqIO.parse(open(input1), 'fasta'):
        id, seq = fasta.description, str(fasta.seq)
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        size = len(seq.replace('N', '').replace('-', ''))
        min_size = genome_size - int(genome_size * max_gaps/100)
        if size > min_size:
            # print(size, min_size)
            # print(id)
            fasta_headers.append(id)
        else:
            print('size: ' + str(size) + ' bp ' + ' - ' + id + ' contains more than ' + str(max_gaps) + '% of Ns. Skipping...')

    print('\n### Loading matrices...')
    # open metadata file
    dfM = pd.read_csv(input2, encoding='utf-8', sep=separator, dtype=str)
    # print(dfM)

    # open genome sampling matrix
    dfS = pd.read_csv(input3, encoding='utf-8', sep=separator, dtype=str)

    # remove genomes with incomplete dates
    dfM = dfM[dfM[date_col].apply(lambda x: len(x.split('-')) == 3)]  # accept only full dates
    dfM = dfM[dfM[date_col].apply(lambda x: 'X' not in x)]  # exclude -XX-XX missing dates


    # filter by date
    today = time.strftime('%Y-%m-%d', time.gmtime())
    dfM[filter_col] = pd.to_datetime(dfM[filter_col]) # converting to datetime format
    if start_date == None:
        start_date = dfM[filter_col].min()
    if end_date == None:
        end_date = today
    # converting back to string
    dfM[filter_col] = dfM[filter_col].apply(lambda x: x.strftime('%Y-%m-%d'))
    # print(dfM)


    # filter genomes based on sampling date
    def filter_bydate(df, date):
        df[date] = pd.to_datetime(df[date]) # converting to datetime format
        # print(df[date])
        mask = (df[date] > start_date) & (df[date] <= end_date) # mask any lines with dates outside the start/end dates
        df = df.loc[mask] # apply mask
        # print(df[date])
        return df
    dfM = filter_bydate(dfM, filter_col)
    # print(dfM)

    # list of sequences to be ignored in all instances
    remove_sequences = []
    for id in open(remove, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            remove_sequences.append(id)

    # check if available sequences have metadata
    meta_seqs = dfM['strain'].to_list()
    intersection = [] # sequences that also have metadata
    for strain in fasta_headers:
        if strain in meta_seqs:
            intersection.append(strain)
            # print('+ ', strain)
        else:
            # print('- ', strain)
            remove_sequences.append(strain)
    # print(intersection)

    # list of sequences to be kept in all instances
    keep_sequences = []
    for id in open(keep, "r").readlines():
        if id[0] not in ["#", "\n"]:
            id = id.strip()
            if id in meta_seqs:
                keep_sequences.append(id)
            else:
                remove_sequences.append(id)

    # keep or remove specific sequences
    dfM = dfM[dfM['strain'].isin(intersection)]
    dfM = dfM[~dfM['strain'].isin(remove_sequences)] # remove bad quality sequences

    # print(keep_sequences)
    # print(remove_sequences)


    # create epiweek column
    def get_epiweeks(date):
        date = pd.to_datetime(date)
        epiweek = str(Week.fromdate(date, system="cdc")) # get epiweeks
        epiweek = epiweek[:4] + '_' + 'EW' + epiweek[-2:]
        return epiweek
    dfM['epiweek'] = dfM[date_col].apply(lambda x: get_epiweeks(x))
    # dfM.set_index('strain', inplace=True)
    # print(dfM)

    # empty matrix dataframe
    columns = sorted(dfM['epiweek'].unique().tolist())
    # print(columns)
    # for i in dfM[geo_level].unique().tolist():
    #     print(i)
    rows = sorted(dfM[geo_level].astype(str).unique().tolist())

    dfG = pd.DataFrame(index=rows, columns=columns)
    dfG.index.name = geo_level
    for column in columns:
        for row in rows:
            dfG.at[row, column] = []

    # print(dfG)
    print(keep_sequences)

    # add pre-selected genomes to matrix
    for genome in keep_sequences:
        if genome in dfM['strain'].to_list():
            # print(genome)
            metadata = dfM.loc[dfM['strain'] == genome]
            location = metadata[geo_level].values[0]
            epiweek = metadata['epiweek'].values[0]
            # print(location, epiweek)
            if genome not in dfG.loc[location, epiweek]:
                dfG.at[location, epiweek] += [genome]
    # print(dfG)
    # print(seed)

    print('\n### Starting sampling process...\n')
    # sampling process
    random.seed(seed) # pseudo-random sampling seed
    glevel = dfM.groupby(geo_level)
    for name, dfLevel in glevel:
        if name in dfS[geo_level].to_list():
            gEpiweek = dfLevel.groupby('epiweek') # geolevel-specific dataframe
            # total_genomes = dfL[geo_level].count() #
            for epiweek, dfEpiweek in gEpiweek:
                available_samples = dfEpiweek['epiweek'].count()  # genomes in bin
                try:
                    target_sampling = int(dfS.loc[dfS[geo_level] == name, epiweek])
                except:
                    target_sampling = available_samples
                # print('')
                # bias = available_samples / target_sampling
                # print(name, epiweek, '-', available_samples, target_sampling, bias)

                existing = dfG.loc[name, epiweek]
                # print(existing)

                if target_sampling >= available_samples:  # if requested sample number is higher than available genomes, get all
                    sampled = [sample for sample in dfEpiweek['strain'].to_list() if sample not in existing]
                    # print(sampled, len(sampled))
                else:
                    pool = [sample for sample in dfEpiweek['strain'].to_list() if sample not in existing]
                    sampled = random.sample(pool, target_sampling - len(existing))
                    # print(sampled, len(sampled))

                dfG.at[name, epiweek] += sampled # add selected samples to dataframe
        # print(dfG.values.tolist())


    # export output
    selected_samples = []
    report = {}
    total_genomes = 0
    # print(dfG)
    for idx, row in dfG.stack().iteritems():
        place = idx[0]
        week = idx[1]
        available = str(len(row))
        try:
            target = str(dfS.loc[dfS[geo_level] == idx[0], idx[1]].values[0])
        except:
            target = available
        print(place + ' on ' + week + ': ' + 'requested = ' + target + '; ' + 'sampled = ' + available)
        # print(idx, dfS.loc[dfS[geo_level] == idx[0], idx[1]].values[0], row)

        name = idx[0]
        if len(row) > 0:
            total_genomes += len(row)
            if name not in report:
                report[name] = 0
            for sample in row:
                selected_samples.append(sample)
            report[name] += len(row)

    # print(report)

    # export fasta file
    print('\n### Exporting sequences and metadata...\n')
    outfile1 = open(output1, 'w')
    c = 1
    found = []
    for fasta in SeqIO.parse(open(input1), 'fasta'):
        id, seq = fasta.description, fasta.seq
        id = id.replace('hCoV-19/', '').split('|')[0].replace(' ', '')
        if id in selected_samples and id not in found:
            print(str(c) + '. ' + id)
            entry = ">" + id + "\n" + str(seq.upper()) + "\n"
            outfile1.write(entry)
            found.append(id)
            c += 1

    # export metadata
    dfM = dfM[dfM['strain'].isin(selected_samples)]
    dfM.to_csv(output2, sep='\t', index=False)


    # export report
    outfile3 = open(output3, 'w')
    outfile3.write('# Seed for pseudo-random sampling: ' + str(seed) + '\n\n')
    outfile3.write('# A total of ' + str(total_genomes) + ' sequences selected from ' + str(len(report)) + ' locations\n\n')
    for loc, count in report.items():
        outfile3.write(str(count) + '\t' + loc + '\n')
    outfile3.write('\n\n# A total of ' + str(len(remove_sequences)) + ' were removed due to lack of metadata, or as listed in remove.txt\n\n')
    for sample in remove_sequences:
        outfile3.write(sample + '\n')
    outfile3.write('\n\n# A total of ' + str(len(keep_sequences)) + ' samples were forcibly added as listed in keep.txt\n\n')
    for sample in keep_sequences:
        outfile3.write(sample + '\n')

    print('\nTotal sampled genomes: ' + str(total_genomes) + '\n')
