#!/usr/bin/python
# -*- coding: utf-8 -*-

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python: 3
#
#  get_daily_matrix_global.py.py -> This code converts Johns Hopkins (US)
#                                   dashboard raw data in CSV into a TSV,
#                                   with reformatted dates and columns.
#
#
# Release date: 2020-09-16
# Last update: 2020-12-08
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import os
import pandas as pd
import argparse
import time


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter newly sequenced genomes not yet added in an existing FASTA dataset of sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=False, help="Johns Hopkins, US daily case counts file")
    parser.add_argument("--download", required=True, nargs=1, type=str, default='yes',
                        choices=['yes', 'no'], help="Download US case series from Johns Hopkins University?")
    parser.add_argument("--start-date", required=False, type=str,  help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str,  help="End date in YYYY-MM-DD format")
    args = parser.parse_args()

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_bubble/nextstrain/run12_20201008_batch15/sampling_prop/data/'
    # input = 'time_series_covid19_confirmed_US.csv'
    # start_date = '2019-12-24' # start date of period of interest
    # end_date = '2020-09-30' # end date of period of interest
    # separator = ','

    download = args.download[0]

    path = os.getcwd()
    if download == 'yes':
        if 'time_series_covid19_usa.csv' not in os.listdir(path):
            url = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv'
            os.system('wget %s' % url)
            os.system('mv %s time_series_covid19_usa.csv' % url.split('/')[-1])
        df = pd.read_csv('time_series_covid19_usa.csv', encoding='utf-8', sep=',', dtype='str')
    else:
        input = args.input
        df = pd.read_csv(input, encoding='utf-8', sep=',', dtype='str')

    start_date = args.start_date
    end_date = args.end_date

    remove = ['Diamond Princess', 'Grand Princess']
    fix_names = {}

    # list of date columns
    date_columns = [column for column in df.columns.to_list() if column[0].isdecimal()]

    # dates = [field for field in df.columns if field[0].isdigit]
    # df = df.drop(columns=['Lat', 'Long']) # drop unwanted columns
    df = df[['Country_Region', 'Province_State'] + date_columns]
    df.fillna('', inplace=True)
    df.insert(2, 'code', '')
    
    # rename column
    df = df.rename(columns={'Country_Region': 'country'})
    df = df.rename(columns={'Province_State': 'state'})
    
    # print(df[['state', 'Lat', 'Long']][df['Lat'] == str(0)])
    
    # drop unwanted rows
    df = df[~df['country'].isin(remove)]
    df = df[~df['state'].isin(remove)]
    
    # rename unusual country and territory names
    df['country'].replace(fix_names, inplace=True)
    df['state'].replace(fix_names, inplace=True)


    today = time.strftime('%Y-%m-%d', time.gmtime())
    if start_date == None:
        start_date = pd.to_datetime([col for col in df.columns.to_list() if col[0].isdecimal()]).min()
    if end_date == None:
        end_date = today

    
    # rename column names and drop columns out of date range
    for column in df.columns.to_list():
        if column[0].isdecimal():
            date = pd.to_datetime(column)
            if date >= pd.to_datetime(start_date) and date <= pd.to_datetime(end_date):
                new_column = date.strftime('%Y-%m-%d')
                df = df.rename(columns={column: new_column})
                df[new_column] = df[new_column].astype(int)
            else:
                df = df.drop(columns=[column])
    
    # print(df.columns)
    
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
        'Northern Mariana Islands':'MP',
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
    
    # add iso code
    df['code'] = df['state'].apply(lambda x: us_state_abbrev[x] if x in us_state_abbrev else 'XX')

    
    # group by country name, summing up values
    df = df.groupby(['code', 'state'], as_index=False).sum()
    

    # chronologically inverted list of dates
    date_columns = [column for column in df.columns.to_list() if column[0].isdecimal()][::-1]
    
    for num, col in enumerate(date_columns):
        # print(num, col)
        if num < len(date_columns)-1: # go up to first day
            for idx, row_value in df[col].iteritems():
                state = df.loc[idx, 'code']
                daily_count = int(row_value) - int(df[date_columns[num+1]][idx])
                if daily_count < 0:
                    daily_count = 0
                df.loc[idx, date_columns[num]] = daily_count
                print(state, date_columns[num], daily_count)#int(row_value), '-', int(df[date_columns[num+1]][idx]), '=', daily_count)
    
    # print(df.head)
    
    # save processed metadata
    df.to_csv('time_series_covid19_usa_reformatted.tsv', sep='\t', index=False)
    print('\nOutput successfully exported: ' + 'time_series_covid19_usa_reformatted.tsv\n')
