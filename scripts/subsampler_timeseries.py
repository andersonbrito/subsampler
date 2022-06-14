import pandas as pd
from Bio import SeqIO
from epiweeks import Week
import random
import time
import argparse
import pycountry_convert as pyCountry
import pycountry
import os
from os import listdir

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=False, help="FASTA file with genomes named as in metadata")
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--genome-matrix", required=True, help="TSV file showing corrected genome counts per unit of time")
    parser.add_argument("--max-missing", required=False, type=int, default=99, help="Maximum percentage of Ns or gaps (int = 1-100)")
    parser.add_argument("--refgenome-size", required=False, type=int, default=1, help="Reference genome size")
    parser.add_argument("--keep", required=False, help="List of samples to keep, in all instances")
    parser.add_argument("--remove", required=False, help="List of samples to remove, in all instances")
    parser.add_argument("--filter-file", required=False, help="TSV file listing columns/values of samples to be batch included or excluded")
    parser.add_argument("--seed", required=False, type=int, help="Seed number for pseudorandom sampling of genomes")
    parser.add_argument("--index-column", required=True, help="Metadata column with unique genome identifiers (genome names, accession codes, etc")
    parser.add_argument("--geo-column", required=True, help="Metadata column with the target geographic information (country, division, etc)")
    parser.add_argument("--date-column", required=True, type=str, help="Metadata column containing the collection dates")
    parser.add_argument("--time-unit", required=True, nargs=1, type=str, default='week', choices=['week', 'month', 'year', 'full'], help="Time unit for conversion")
    parser.add_argument("--weekasdate",required=False, nargs=1, type=str, default='no', choices=['start', 'end', 'no'], help="When representing weeks as date, which day of the week should be used?")
    parser.add_argument("--start-date", required=False, type=str, help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str, help="End date in YYYY-MM-DD format")
    parser.add_argument("--sampled-sequences", required=True, help="Sampled genomes")
    parser.add_argument("--sampled-metadata", required=True, help="Sampled metadata")
    parser.add_argument("--report", required=True, help="List of statistics related to the sampling scheme")
    args = parser.parse_args()

    metadata_file = args.metadata
    sampling_file = args.genome_matrix
    fasta_file = args.sequences
    genome_size = args.refgenome_size
    max_gaps = args.max_missing
    keep = args.keep
    remove = args.remove
    filter_file = args.filter_file
    seed = args.seed
    id_col = args.index_column
    geo_level = args.geo_column
    date_col = args.date_column
    unit = args.time_unit[0]
    weekasdate = args.weekasdate[0]
    start_date = args.start_date
    end_date = args.end_date
    outfile_sequences = args.sampled_sequences
    outfile_metadata = args.sampled_metadata
    outfile_report = args.report

    # path = '/Users/anderson/Desktop/subsampler_issues/20220611_update/'
    # os.chdir(path)
    # metadata_file = 'data/metadata.tsv'
    # sampling_file = 'outputs/matrix_genomes_unit_corrected.tsv'
    # fasta_file = '' #path + 'data/sequences.fasta'
    # filter_file = 'config/filters.tsv'
    # # include_file = path + 'config/strict_inclusion.tsv'
    # # drop_file = path + 'config/batch_removal.tsv'
    # keep = 'config/keep.txt'
    # remove = 'config/remove.txt'
    # id_col = 'gisaid_epi_isl'
    # geo_level = 'country_exposure'
    # seed = 2007
    # # seed = None
    # date_col = 'date'
    # # filter_col = 'date'
    # genome_size = 29930
    # max_gaps = 30
    # unit = 'month'
    # weekasdate = 'no'
    # start_date = '2020-03-01'
    # end_date = '2021-01-31'
    # # start_date = None
    # # end_date = None
    # outfile_sequences = 'sequences_corrected.txt'
    # outfile_metadata = 'metadata_corrected.tsv'
    # outfile_report = 'report.txt'

    if seed == None:
        seed = random.random()

    # if filter_col == None:
    #     filter_col = date_col

    # pd.set_option('display.max_columns', 500)

    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    print('\n### Loading matrices...')
    # open metadata file
    dfM = load_table(metadata_file)
    dfM.fillna('', inplace=True)
    # print(dfM)


    # get sequence headers
    fasta_headers = []
    if fasta_file not in ['', None]:
        print('\n### Loading sequences...')
        for fasta in SeqIO.parse(open(fasta_file), 'fasta'):
            id, seq = fasta.description, str(fasta.seq)
            id = id.split('|')[0].replace(' ', '')
            size = len(seq.replace('N', '').replace('-', ''))
            min_size = genome_size - int(genome_size * max_gaps / 100)
            if size > min_size:
                fasta_headers.append(id)
            else:
                print('size: ' + str(size) + ' bp ' + ' - ' + id + ' contains more than ' + str(
                    max_gaps) + '% of Ns. Skipping...')
    else:
        fasta_headers = list(set(dfM[id_col].tolist()))


    # filter rows
    def filter_df(df, criteria):
        print('\n### Filtering rows...')
        new_df = pd.DataFrame()
        include = {}
        for filter_value in criteria.split(','):
            filter_value = filter_value.strip()
            if not filter_value.startswith('~'):
                col, val = filter_value.split(':')[0], filter_value.split(':')[1]
                if val == '\'\'':
                    val = ''
                if col not in include:
                    include[col] = [val]
                else:
                    include[col].append(val)
        # print('Include:', include)
        for filter_col, filter_val in include.items():
            print('\t- Including only rows with \'' + filter_col + '\' = \'' + ', '.join(filter_val) + '\'')
            # print(new_df.size)
            if new_df.empty:
                df_filtered = df[df[filter_col].isin(filter_val)]
                new_df = new_df.append(df_filtered)
            else:
                new_df = new_df[new_df[filter_col].isin(filter_val)]
            # print(new_df)#.head())

        exclude = {}
        for filter_value in criteria.split(','):
            filter_value = filter_value.strip()
            if filter_value.startswith('~'):
                # print('\t- Excluding all rows with \'' + col + '\' = \'' + val + '\'')
                filter_value = filter_value[1:]
                col, val = filter_value.split(':')[0], filter_value.split(':')[1]
                if val == '\'\'':
                    val = ''
                if col not in exclude:
                    exclude[col] = [val]
                else:
                    exclude[col].append(val)
        # print('Exclude:', exclude)
        for filter_col, filter_val in exclude.items():
            print('\t- Excluding all rows with \'' + filter_col + '\' = \'' + ', '.join(filter_val) + '\'')
            if new_df.empty:
                df = df[~df[filter_col].isin(filter_val)]
                new_df = new_df.append(df)
            else:
                new_df = new_df[~new_df[filter_col].isin(filter_val)]
            # print(new_df)#.head())
        return new_df


    # filtering criteria
    if filter_file not in ['', None]:
        dfC = load_table(filter_file)
        dfC['action'] = dfC['action'].apply(lambda x: '~' if x == 'exclude' else '')  # exclude -XX-XX missing dates
        dfC['filter'] = dfC['action'].astype(str) + dfC['column'].astype(str) + ':' + dfC['value'].astype(str)
        filters = ', '.join(dfC['filter'].tolist())
        dfM = filter_df(dfM, filters)


    # open genome sampling matrix
    dfS = load_table(sampling_file)

    print('\n### Removing genomes with incomplete dates')
    # remove genomes with incomplete dates
    dfM = dfM[dfM[date_col].apply(lambda x: len(x.split('-')) == 3)]  # accept only full dates
    dfM = dfM[dfM[date_col].apply(lambda x: 'X' not in x)]  # exclude -XX-XX missing dates


    # filter by date
    today = time.strftime('%Y-%m-%d', time.gmtime())
    dfM[date_col] = pd.to_datetime(dfM[date_col])  # converting to datetime format
    if start_date == None:
        start_date = dfM[date_col].min()
    if end_date == None:
        end_date = today

    # converting dates back to string format
    dfM[date_col] = dfM[date_col].apply(lambda x: x.strftime('%Y-%m-%d'))


    print('\n### Filtering genomes by date')
    # filter genomes based on sampling date
    def filter_bydate(df, date):
        df[date] = pd.to_datetime(df[date])  # converting to datetime format
        mask = (df[date] > start_date) & (df[date] <= end_date)  # mask any lines with dates outside the start/end dates
        df = df.loc[mask]  # apply mask
        return df

    dfM = filter_bydate(dfM, date_col)
    # print(dfM)

    print('\n### Removing genomes tagged for removal')
    # list of sequences to be ignored in all instances
    remove_sequences = []
    if remove not in ['', None]:
        for id in open(remove, "r").readlines():
            if id[0] not in ["#", "\n"]:
                id = id.strip()
                remove_sequences.append(id)

    # print('\t- Checking if genomes have metadata, and removing if negative')

    # check if available sequences have metadata
    meta_seqs = dfM[id_col].to_list()
    intersection = set(fasta_headers).intersection(meta_seqs)

    def Diff(li1, li2):
        return (list(list(set(li1) - set(li2)) + list(set(li2) - set(li1))))


    remove_sequences = remove_sequences + Diff(intersection, meta_seqs)

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
    dfM = dfM[dfM[id_col].isin(intersection)]  # include only sequences with metadata
    dfM = dfM[~dfM[id_col].isin(remove_sequences)]  # remove bad quality sequences


    ### FIX OR ADD NEW COLUMNS IN THE METADATA

    # converting dates back to string format
    dfM[date_col] = dfM[date_col].apply(lambda x: x.strftime('%Y-%m-%d'))

    # create time unit column
    # time_cols = []
    def get_newunit(value):
        if value[0].isdecimal():
            date = pd.to_datetime(value)
            if unit == 'week':
                epiweek = str(Week.fromdate(date, system="cdc")) # get epiweeks
                year, week = epiweek[:4], epiweek[-2:]
                if weekasdate in ['start', 'end']:
                    if weekasdate == 'start':
                        epiweek = str(Week(int(year), int(week)).startdate())
                    else:
                        epiweek = str(Week(int(year), int(week)).enddate())
                else:
                    epiweek = year + '_' + 'EW' + week
                # if epiweek not in time_cols:
                #     time_cols.append(epiweek)
                return epiweek
            elif unit == 'month':
                year_month = date.strftime("%Y-%m")
                # if year_month not in time_cols:
                    # time_cols.append(year_month)
                return year_month
            elif unit == 'year':
                year = date.strftime("%Y")
                # if year not in time_cols:
                    # time_cols.append(year)
                return year
            elif unit == 'full':
                return 'total'
        else:
            if unit == 'full':
                return 'total'
            else:
                return value

    dfM['time_unit'] = dfM[date_col].apply(lambda x: get_newunit(x))



    # fix place of origin when disagreements between 'place' and 'place_exposure' exist
    if 'exposure' in geo_level:
        geo_columns = ['region', 'country', 'division']
        for level in geo_columns:
            exposure_column = level + '_exposure'
            for idx, row in dfM.iterrows():
                if dfM.loc[idx, exposure_column].lower() in ['', 'unknown']:
                    dfM.loc[idx, exposure_column] = dfM.loc[idx, level]


    # get ISO alpha3 country codes
    codes = {}
    def get_iso(country):
        global codes
        if country not in codes.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                codes[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    codes[country] = isoCode
                except:
                    codes[country] = ''
        return codes[country]


    us_state_abbrev = {
        'Alabama': 'AL',
        'Alaska': 'AK',
        'American Samoa': 'AS',
        'Arizona': 'AZ',
        'Arkansas': 'AR',
        'California': 'CA',
        'Colorado': 'CO',
        'Connecticut': 'CT',
        'Delaware': 'DE',
        'District of Columbia': 'DC',
        'Washington DC': 'DC',
        'Florida': 'FL',
        'Georgia': 'GA',
        'Guam': 'GU',
        'Hawaii': 'HI',
        'Idaho': 'ID',
        'Illinois': 'IL',
        'Indiana': 'IN',
        'Iowa': 'IA',
        'Kansas': 'KS',
        'Kentucky': 'KY',
        'Louisiana': 'LA',
        'Maine': 'ME',
        'Maryland': 'MD',
        'Massachusetts': 'MA',
        'Michigan': 'MI',
        'Minnesota': 'MN',
        'Mississippi': 'MS',
        'Missouri': 'MO',
        'Montana': 'MT',
        'Nebraska': 'NE',
        'Nevada': 'NV',
        'New Hampshire': 'NH',
        'New Jersey': 'NJ',
        'New Mexico': 'NM',
        'New York': 'NY',
        'North Carolina': 'NC',
        'North Dakota': 'ND',
        'Northern Mariana Islands': 'MP',
        'Ohio': 'OH',
        'Oklahoma': 'OK',
        'Oregon': 'OR',
        'Pennsylvania': 'PA',
        'Puerto Rico': 'PR',
        'Rhode Island': 'RI',
        'South Carolina': 'SC',
        'South Dakota': 'SD',
        'Tennessee': 'TN',
        'Texas': 'TX',
        'Utah': 'UT',
        'Vermont': 'VT',
        'Virgin Islands': 'VI',
        'Virginia': 'VA',
        'Washington': 'WA',
        'West Virginia': 'WV',
        'Wisconsin': 'WI',
        'Wyoming': 'WY'
    }

    # add state code
    if 'code' not in dfM.columns.to_list():
        dfM.insert(1, 'code', '')
        if 'division' in geo_level:
            dfM['code'] = dfM[geo_level].apply(lambda x: us_state_abbrev[x] if x in us_state_abbrev else '')
        elif 'country' in geo_level:
            dfM['code'] = dfM[geo_level].apply(lambda x: get_iso(x))
        else:
            dfM['code'] = dfM[geo_level]

    # set geo_level as code
    geo_level = 'code'

    # empty matrix dataframe
    columns = sorted(dfM['time_unit'].unique().tolist())
    rows = sorted(dfM[geo_level].astype(str).unique().tolist())

    dfG = pd.DataFrame(index=rows, columns=columns)
    dfG.index.name = geo_level
    for column in columns:
        for row in rows:
            dfG.at[row, column] = []


    # add pre-selected genomes to matrix
    for genome in keep_sequences:
        if genome in dfM[id_col].to_list():
            # print(genome)
            metadata = dfM.loc[dfM[id_col] == genome]
            location = metadata[geo_level].values[0]
            time_unit = metadata['time_unit'].values[0]
            # print(location, time_unit)
            if genome not in dfG.loc[location, time_unit]:
                dfG.at[location, time_unit] += [genome]
    # print(dfG)
    # print(seed)


    print('\n### Starting sampling process...\n')
    # sampling process
    random.seed(seed)  # pseudo-random sampling seed
    glevel = dfM.groupby(geo_level)
    for name, dfLevel in glevel:
        if name in dfS[geo_level].to_list():
            gUnit = dfLevel.groupby('time_unit')  # geolevel-specific dataframe
            for time_unit, dfUnit in gUnit:
                available_samples = dfUnit['time_unit'].count()  # genomes in bin
                try:
                    target_sampling = int(dfS.loc[dfS[geo_level] == name, time_unit])
                except:
                    target_sampling = 1  # available_samples # take this number of genomes when not epidata is available
                # print('')
                # print(name, time_unit, '-', available_samples, target_sampling, bias)

                existing = dfG.loc[name, time_unit]  # pre-selected sequences, if any was listed in keep.txt

                if target_sampling >= available_samples:  # if requested sample number is higher than available genomes, get all
                    sampled = [sample for sample in dfUnit[id_col].to_list() if sample not in existing]
                    # print(sampled, len(sampled))
                elif target_sampling == 1:
                    pool = [sample for sample in dfUnit[id_col].to_list() if sample not in existing]
                    sampled = random.sample(pool, 1)
                else:
                    pool = [sample for sample in dfUnit[id_col].to_list() if sample not in existing]
                    if target_sampling < len(existing):
                        target_sampling = len(existing)
                    sampled = random.sample(pool, target_sampling - len(existing))

                dfG.at[name, time_unit] += sampled  # add selected samples to dataframe

    # export output
    selected_samples = []
    report = {}
    total_genomes = 0

    for idx, row in dfG.stack().iteritems():
        place = idx[0]
        time_period = idx[1]
        available = str(len(row))
        try:
            target = str(dfS.loc[dfS[geo_level] == idx[0], idx[1]].values[0])
        except:
            target = available
        print('\t- ' + place + ' on ' + time_period + ': ' + 'requested = ' + target + '; ' + 'sampled = ' + available)

        name = idx[0]
        if len(row) > 0:
            total_genomes += len(row)
            if name not in report:
                report[name] = 0
            for sample in row:
                selected_samples.append(sample)
            report[name] += len(row)

    # export fasta file
    print('\n### Exporting sequence list and metadata...\n')
    outfile1 = open(outfile_sequences, 'w')
    c = 1
    found = []
    for id in selected_samples:
        if id not in found:
            print('\t- ' + str(c) + '. ' + id)
            outfile1.write(id + '\n')
            found.append(id)
            c += 1

    # export metadata
    dfM = dfM[dfM[id_col].isin(selected_samples)]
    dfM = dfM.sort_values(by=geo_level)
    dfM.to_csv(outfile_metadata, sep='\t', index=False)

    # export report
    outfile3 = open(outfile_report, 'w')
    outfile3.write('# Seed for pseudo-random sampling: ' + str(seed) + '\n\n')
    outfile3.write(
        '# A total of ' + str(total_genomes) + ' sequences selected from ' + str(len(report)) + ' locations\n\n')
    for loc, count in report.items():
        outfile3.write(str(count) + '\t' + loc + '\n')
    outfile3.write('\n\n# A total of ' + str(
        len(remove_sequences)) + ' were removed due to lack of metadata, or as listed in remove.txt\n\n')

    for sample in remove_sequences:
        outfile3.write(sample + '\n')
    outfile3.write('\n\n# A total of ' + str(len(keep_sequences)) + ' samples were forcibly added as listed in keep.txt\n\n')

    for sample in keep_sequences:
        outfile3.write(sample + '\n')

    print('\nTotal sampled genomes: ' + str(total_genomes) + '\n')

