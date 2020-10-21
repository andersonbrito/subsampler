#!/usr/bin/python
# -*- coding: utf-8 -*-

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python: 3
#
#   .py -> This code reads DBF files from DataSUS (Brazil's
#                       Ministry of Health), group state level data as
#                       national level case counts, and export them to
#                       a CSV file.
#
# Release date: 14/Feb/2020
# Last update: 14/Feb/2020
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import os
import pandas as pd
import pycountry_convert as pyCountry
import pycountry
import argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter newly sequenced genomes not yet added in an existing FASTA dataset of sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="Johns Hopkins, global daily case counts file")
    parser.add_argument("--start-date", required=False, type=str,  help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end-date", required=False, type=str,  help="End date in YYYY-MM-DD format")
    args = parser.parse_args()

    input = args.input
    start_date = args.start_date
    end_date = args.end_date
    separator = ','


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_bubble/nextstrain/run12_20201008_batch15/sampling_prop/data/'
    # input = 'time_series_covid19_confirmed_global.csv'
    # start_date = '2019-12-24' # start date of period of interest
    # end_date = '2020-09-30' # end date of period of interest

    remove = ['Diamond Princess', 'MS Zaandam', 'Channel Islands']
    fix_names = {'Congo (Brazzaville)': 'Congo', 'Congo (Kinshasa)': 'Democratic Republic of the Congo',
                 'St Martin': 'Saint Martin', 'Korea, South': 'South Korea', 'Taiwan*': 'Taiwan',
                 'West Bank and Gaza': 'Palestine', 'Burma': 'Myanmar', '"Bonaire, Sint Eustatius and Saba"': 'Bonaire'}


    df = pd.read_csv(input, encoding='utf-8', sep=separator, dtype=str)

    # list of date columns
    date_columns = [column for column in df.columns.to_list() if column[0].isdecimal()]
    df = df[['Country/Region', 'Province/State'] + date_columns]
    df.fillna('', inplace=True)
    df.insert(2, 'iso', '')


    # drop unwanted rows
    df = df[~df['Country/Region'].isin(remove)]
    df = df[~df['Province/State'].isin(remove)]

    print(df)

    # rename unusual country and territory names
    df['Country/Region'].replace(fix_names, inplace=True)
    df['Province/State'].replace(fix_names, inplace=True)


    # rename column names and drop columns out of date range
    today = time.strftime('%Y-%m-%d', time.gmtime())
    if start_date == None:
        start_date = pd.to_datetime([col for col in df.columns.to_list() if col[0].isdecimal()]).min()
    if end_date == None:
        end_date = today

    for column in df.columns.to_list():
        if column[0].isdecimal():
            date = pd.to_datetime(column)
            if date >= pd.to_datetime(start_date) and date <= pd.to_datetime(end_date):
                new_column = date.strftime('%Y-%m-%d')
                df = df.rename(columns={column: new_column})
                df[new_column] = df[new_column].astype(int)
            else:
                df = df.drop(columns=[column])


    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # add iso code
    df['iso'] = df['Country/Region'].map(get_iso)

    # country dependencies
    has_dependencies = ['United Kingdom', 'France', 'Denmark', 'Netherlands']
    is_autonomous = ['Hong Kong']#, 'Curacao', 'Sint Maarten', 'Bonaire']
    for idx, row in df.iterrows():
        country = df.loc[idx, 'Country/Region']
        province = df.loc[idx, 'Province/State']
        iso = df.loc[idx, 'iso']
        if province in is_autonomous:
            df.loc[idx, 'Country/Region'] = province
            df.loc[idx, 'iso'] = get_iso(province)

        if province not in '' and country in has_dependencies:
            # print(country, province)
            if get_iso(province) == get_iso(country):
                continue
            else:
                df.loc[idx, 'Country/Region'] = province
                df.loc[idx, 'iso'] = get_iso(province)


    df = df.drop(columns=['Province/State']) # drop unwanted columns
    df = df.rename(columns={'Country/Region': 'country'})


    # group by country name, summing up values
    df = df.groupby(['iso', 'country'], as_index=False).sum()


    date_columns = [column for column in df.columns.to_list() if column[0].isdecimal()][::-1]

    # print(date_columns)
    # print(date_columns[::-1])

    # convert cumulative counts into daily counts
    for num, col in enumerate(date_columns):
        # print(num, col)
        if num < len(date_columns)-1:
            # print(num, len(date_columns), date_columns[len(date_columns)-1])
            for idx, row_value in df[col].iteritems():
                country = df.loc[idx, 'iso']
                daily_count = int(row_value) - int(df[date_columns[num+1]][idx])
                if daily_count < 0:
                    daily_count = 0
                df.loc[idx, date_columns[num]] = daily_count
                print(country, date_columns[num], daily_count)#, int(row_value), int(df[date_columns[num+1]][idx]))


    # save processed metadata
    df.to_csv(input.split('.')[0] + '_reformatted.tsv', sep='\t', index=False)
    print('\nOutput successfully exported: ' + input.split('.')[0] + '_reformatted.tsv\n')