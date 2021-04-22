import pandas as pd
import time
import argparse
import pycountry_convert as pyCountry
import pycountry

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--index-column", required=True, help="Column with unique geographic information")
    parser.add_argument("--extra-columns", required=False, nargs='+', type=str,
                        help="extra columns with geographic info to export")
    parser.add_argument("--date-column", required=True, type=str, help="Column containing the date information")
    parser.add_argument("--start-date", required=False, type=str, help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str, help="End date in YYYY-MM-DD format")
    parser.add_argument("--output", required=True, help="Genome matrix")
    args = parser.parse_args()

    metadata = args.metadata
    geo_col = args.index_column
    extra_cols = args.extra_columns
    date_col = args.date_column
    group_by = ['code', date_col]
    start_date = args.start_date
    end_date = args.end_date
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/run15_20210422_samprop/'
    # metadata = path + 'data/metadata_nextstrain.tsv'
    # output = path + 'matrix_genomes_daily.tsv'

    # geo_col = 'division_exposure'
    # date_col = 'date'
    # extra_cols = ['country_exposure']
    # group_by = ['code', date_col]
    # start_date = '2019-12-01'
    # end_date = '2020-07-22'
    # start_date = None
    # end_date = None

    pd.set_option('display.max_columns', 500)

    # input genome and case counts per epiweek
    df = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    df.fillna('', inplace=True)

    # fix exposure
    geolevels = ['region', 'country', 'division']
    print('\n * Loading genome metadata\n')
    for level in geolevels:
        exposure_column = level + '_exposure'
        if exposure_column == geo_col:
            for idx, row in df.iterrows():
                if df.loc[idx, exposure_column].lower() in ['', 'unknown']:
                    df.loc[idx, exposure_column] = df.loc[idx, level]

    # get ISO alpha3 country codes
    codes = {'Rest of the US': 'RES', 'NewYork/NewJersey': 'NYJ'}
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
    print('\n * Converting ' + geo_col + ' into codes (acronyms)\n')
    if 'code' not in df.columns.to_list():
        df.insert(1, 'code', '')
        if 'division' in geo_col:
            df['code'] = df[geo_col].apply(lambda x: us_state_abbrev[x] if x in us_state_abbrev else '')
        elif 'country' in geo_col:
            df['code'] = df[geo_col].apply(lambda x: get_iso(x))
        else:
            df['code'] = df[geo_col]

    # remove genomes with incomplete dates
    print('\n * Removing genomes with incomplete dates\n')
    df = df[df[date_col].apply(lambda x: len(x.split('-')) == 3)]  # accept only full dates
    df = df[df[date_col].apply(lambda x: 'X' not in x)]  # exclude -XX-XX missing dates

    # filter by date
    print('\n * Filtering genomes by start and end dates\n')
    today = time.strftime('%Y-%m-%d', time.gmtime())
    df[date_col] = pd.to_datetime(df[date_col])  # converting to datetime format
    if start_date == None:
        start_date = df[date_col].min()
    if end_date == None:
        end_date = today

    mask = (df[date_col] >= start_date) & (
                df[date_col] <= end_date)  # mask any lines with dates outside the start/end dates
    df = df.loc[mask]  # apply mask

    # filter out genomes with missing 'geo_level' name
    df = df[df['code'].apply(lambda x: len(str(x)) > 0)]

    # report
    print('\n* Available genomes\n')
    print('\tOldest collected sampled = ' + df[date_col].min().strftime('%Y-%m-%d'))
    print('\tNewest collected sampled = ' + df[date_col].max().strftime('%Y-%m-%d'))
    print('')

    # convert back to string format
    df[date_col] = df[date_col].apply(lambda x: x.strftime('%Y-%m-%d'))


    # group lines based on date and geolocation, and return genome counts
    df2 = df.groupby(group_by).size().to_frame(name='genome_count').reset_index()
    # print(df2)

    columns = sorted(df[date_col].unique().tolist())
    rows = sorted(df['code'].unique().tolist())

    # empty matrix dataframe
    df3 = pd.DataFrame(index=rows, columns=columns)
    df3 = df3.fillna(0)  # with 0s rather than NaNs

    # give index a name
    df3.index.name = 'code'
    # print(df3)

    # add other columns, if available
    if extra_cols == None:
        extra_cols = []

    for column in extra_cols:
        if column in df.columns.to_list():
            df3.insert(0, column, '')
    df.set_index('code', inplace=True)

    # fill extra columns with their original content
    for idx, row in df3.iterrows():
        for column in extra_cols:
            if column in df.columns.to_list():
                # value = df.loc[idx, column][0]
                value = df.loc[df.index == idx][column].values[0]
                df3.at[idx, column] = value

    # fill matrix with genome counts
    print('\n * Exporting matrix of daily genome counts\n')
    found = []
    for idx, row in df2.iterrows():
        geo = df2.loc[idx, 'code']
        time = df2.loc[idx, date_col]
        count = df2.loc[idx, 'genome_count']
        df3.at[geo, time] = count

    # output processed dataframe
    df3.to_csv(output, sep='\t', index=True)
