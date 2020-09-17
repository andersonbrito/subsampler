import pandas as pd
from epiweeks import Week
import argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="TSV matrix of daily case counts per location")
    parser.add_argument("--start-date", required=False, type=str,  help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str,  help="End date in YYYY-MM-DD format")
    parser.add_argument("--output", required=True, help="TSV matrix of weekly case counts per location")
    args = parser.parse_args()


    input = args.input
    start_date = args.start_date
    end_date = args.end_date
    output = args.output
    separator = '\t'


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_samplingbias/subsampling/run2/'
    # input = path + 'matrix_genomes_days.tsv'
    # output = input.split('.')[0] + '_epiweeks.tsv'
    # separator = '\t'
    # #
    # # start_date = '2019-12-29' # start date of period of interest
    # # end_date = '2020-02-15' # end date of period of interest
    # start_date = None
    # end_date = None


    # open dataframes
    df = pd.read_csv(input, encoding='utf-8', sep=separator, dtype=str)

    # rename column names and drop columns out of date range
    today = time.strftime('%Y-%m-%d', time.gmtime())
    if start_date == None:
        start_date = pd.to_datetime([col for col in df.columns.to_list() if col[0].isdecimal()]).min()
    if end_date == None:
        end_date = today

    nondate_cols = []
    def filterDF(df):
        for column in df.columns.to_list():
            # print(column)
            if column[0].isdecimal():
                date = pd.to_datetime(column)
                if date >= pd.to_datetime(start_date) and date <= pd.to_datetime(end_date):
                    new_column = date.strftime('%Y-%m-%d')
                    df = df.rename(columns={column: new_column})
                    df[new_column] = df[new_column].astype(int)
                else:
                    df = df.drop(columns=[column])
            else:
                if column not in nondate_cols:
                    nondate_cols.append(column)
        return df


    # convert date to epiweek
    ew_cols = []
    def get_epiweeks(value):
        if value[0].isdecimal():
            date = pd.to_datetime(value)
            epiweek = str(Week.fromdate(date, system="cdc")) # get epiweeks
            epiweek = epiweek[:4] + '_' + 'EW' + epiweek[-2:]
            if epiweek not in ew_cols:
                ew_cols.append(epiweek)
            return epiweek
        else:
            return value


    # filter, transpose, convert dates to epiweeks, and re-transpose
    def epiweek_covert(df):
        df = filterDF(df).transpose()
        df['epiweek'] = df.index.map(get_epiweeks) # create new column 'epiweek', mapping 'dates'
        df = df.groupby(['epiweek'], as_index=True).sum() # group dates from same epiweek, sum up counts
        return df.transpose()

    df = epiweek_covert(df)

    df = df[nondate_cols + sorted(ew_cols)]

    # output converted dataframes
    df.to_csv(output, sep='\t', index=False)

