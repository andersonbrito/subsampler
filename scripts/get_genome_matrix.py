import pandas as pd
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--index-column", required=True, help="Column with unique geographic information")
    parser.add_argument("--extra-columns", required=False, nargs='+', type=str,  help="extra columns with geographic info to export")
    parser.add_argument("--date-column", required=True, type=str,  help="Column containing the date information")
    parser.add_argument("--start-date", required=False, type=str,  help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str,  help="End date in YYYY-MM-DD format")
    parser.add_argument("--output", required=True, help="Genome matrix")
    args = parser.parse_args()


    metadata = args.metadata
    geo_col = args.index_column
    extra_cols = args.extra_columns
    # extra_cols = [col.strip() for col in args.extra_columns[0].split()]
    date_col = args.date_column
    group_by = [geo_col, date_col]
    start_date = args.start_date
    end_date = args.end_date
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_samplingbias/subsampling/run3_real/'
    # metadata = path + 'metadata_nextstrain_exposure_iso.tsv'
    # output = path + 'matrix_genomes_daily.tsv'
    #
    # geo_col = 'iso'
    # date_col = 'date'
    # extra_cols = ['region', 'country']
    # group_by = [geo_col, date_col]
    # # start_date = '2019-12-01'
    # # end_date = '2020-02-15'
    # start_date = None
    # end_date = None


    pd.set_option('display.max_columns', 500)

    # input genome and case counts per epiweek
    df = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    df.fillna('', inplace=True)


    # remove genomes with incomplete dates
    df = df[df[date_col].apply(lambda x: len(x.split('-')) == 3)]  # accept only full dates
    df = df[df[date_col].apply(lambda x: 'X' not in x)]  # exclude -XX-XX missing dates

    # filter by date
    today = time.strftime('%Y-%m-%d', time.gmtime())
    df[date_col] = pd.to_datetime(df[date_col]) # converting to datetime format
    if start_date == None:
        start_date = df[date_col].min()
    if end_date == None:
        end_date = today

    # # remove genomes with incomplete dates
    # df[date_col] = df[date_col].apply(lambda x: x.strftime('%Y-%m-%d'))


    mask = (df[date_col] >= start_date) & (df[date_col] <= end_date) # mask any lines with dates outside the start/end dates
    df = df.loc[mask] # apply mask

    # report
    print('\n### Available genomes\n')
    print('Oldest collected sampled = ' + df[date_col].min().strftime('%Y-%m-%d'))
    print('Newest collected sampled = ' + df[date_col].max().strftime('%Y-%m-%d'))
    print('')

    # convert back to string format
    df[date_col] = df[date_col].apply(lambda x: x.strftime('%Y-%m-%d'))

    # filter out genomes with missing 'geo_level' name
    df = df[df[geo_col].apply(lambda x: len(str(x)) > 0)]


    # filter out genomes with incomplete dates
    # df = df[df[date_col].apply(lambda x: len(x.split('-')) == 3)]  # accept only full dates
    # df = df[df[date_col].apply(lambda x: 'X' not in x)]  # exclude -XX-XX missing dates
    # print(df)


    # group lines based on date and geolocation, and return genome counts
    df2 = df.groupby(group_by).size().to_frame(name='genome_count').reset_index()
    # print(df2)

    columns = sorted(df[date_col].unique().tolist())
    rows = sorted(df[geo_col].unique().tolist())

    # empty matrix dataframe
    df3 = pd.DataFrame(index=rows, columns=columns)
    df3 = df3.fillna(0) # with 0s rather than NaNs
    # give index a name
    df3.index.name = geo_col
    # print(df3)


    # add other columns, if available
    for column in extra_cols:
        if column in df.columns.to_list():
            df3.insert(0, column, '')
    df.set_index(geo_col, inplace=True)

    # fill extra columns with their original content
    for idx, row in df3.iterrows():
        for column in extra_cols:
            if column in df.columns.to_list():
                # value = df.loc[idx, column][0]
                value = df.loc[df.index == idx][column].values[0]
                # print(value)
                df3.at[idx, column] = value

    # fill matrix with genome counts
    found = []
    for idx, row in df2.iterrows():
        geo = df2.loc[idx, geo_col]
        time = df2.loc[idx, date_col]
        count = df2.loc[idx, 'genome_count']
        df3.at[geo, time] = count


    # output processed dataframe
    df3.to_csv(output, sep='\t', index=True)
