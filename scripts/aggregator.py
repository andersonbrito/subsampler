import pandas as pd
from epiweeks import Week
import argparse
import time

import platform
# print('Python version:', platform.python_version())
# print('Pandas version:', pd.__version__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Aggregate daily counts as epiweeks, months or year",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="Matrix of daily counts per location")
    parser.add_argument("--unit", required=True, nargs=1, type=str, default='week',
                        choices=['week', 'month', 'year', 'full'], help="Time unit for conversion")
    parser.add_argument("--format",required=False, nargs=1, type=str, default='float',
                        choices=['float', 'integer'], help="What is the format of the data points (float/integer)?")
    parser.add_argument("--weekasdate",required=False, nargs=1, type=str, default='no',
                        choices=['start', 'end', 'no'], help="If representing weeks as date, which day of the week should be used?")
    parser.add_argument("--start-date", required=False, type=str,  help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str,  help="End date in YYYY-MM-DD format")
    parser.add_argument("--output", required=True, help="TSV matrix with aggregated counts")
    args = parser.parse_args()


    input = args.input
    unit = args.unit[0]
    data_format = args.format[0]
    weekasdate = args.weekasdate[0]
    start_date = args.start_date
    end_date = args.end_date
    output = args.output


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/ITpS/projetos_itps/metasurvBR/analyses/bubbles_20211016/'
    # input = path + 'cases_SE35-40_cidades.tsv'
    # unit = 'week'
    # output = input.split('.')[0] + '_' + unit + '.tsv'
    #
    # start_date = '2020-04-01' # start date of period of interest
    # end_date = '2021-05-31' # end date of period of interest
    # start_date = None
    # end_date = None

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

    # Load metadata
    df = load_table(input)

    # rename column names and drop columns out of date range
    today = time.strftime('%Y-%m-%d', time.gmtime())
    if start_date == None:
        start_date = pd.to_datetime([col for col in df.columns.to_list() if col[0].isdecimal()]).min()
    if end_date == None:
        end_date = today

    nondate_cols = []
    def filter_bydate(df):
        for column in df.columns.to_list():
            if column[0].isdecimal():
                date = pd.to_datetime(column)
                if date >= pd.to_datetime(start_date) and date <= pd.to_datetime(end_date):
                    new_column = date.strftime('%Y-%m-%d')
                    df = df.rename(columns={column: new_column})
                    df[new_column] = df[new_column].astype(float)
                    if data_format == 'integer':
                        df[new_column] = df[new_column].astype(int)
                else:
                    df = df.drop(columns=[column])
            else:
                if column not in nondate_cols:
                    nondate_cols.append(column)
        return df

    # convert date
    time_cols = []
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
                if epiweek not in time_cols:
                    time_cols.append(epiweek)
                return epiweek
            elif unit == 'month':
                year_month = date.strftime("%Y-%m")
                if year_month not in time_cols:
                    time_cols.append(year_month)
                return year_month
            elif unit == 'year':
                year = date.strftime("%Y")
                if year not in time_cols:
                    time_cols.append(year)
                return year
            elif unit == 'full':
                return 'total'
        else:
            if unit == 'full':
                return 'total'
            else:
                return value

    # print(df.head())
    # filter, transpose, convert dates to epiweeks, and re-transpose
    def unit_coverter(df):
        df = filter_bydate(df).transpose()
        df['time_variable'] = df.index.map(get_newunit) # create new column 'time_variable', mapping 'dates'
        df = df.groupby(['time_variable'], as_index=True).sum() # group dates from same 'unit', sum up counts
        return df.transpose()

    df = unit_coverter(df)
    df = df[nondate_cols + sorted(time_cols)]

    # output converted dataframes
    df.to_csv(output, sep='\t', index=False)

    print('\nConversion successfully completed.\n')