import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genome-matrix", required=True, help="TSV file showing the original genome counts per epiweek")
    parser.add_argument("--case-matrix", required=True, help="TSV file showing the case counts per epiweek")
    parser.add_argument("--index-column", required=True, help="Column with unique geographic information")
    parser.add_argument("--baseline", required=False, type=float,  help="Baseline sampling proportion")
    parser.add_argument("--output1", required=True, help="TSV file showing genome sampling proportions per epiweek")
    parser.add_argument("--output2", required=True, help="TSV file showing genome sampling bias per epiweek")
    parser.add_argument("--output3", required=True, help="TSV file showing corrected genome counts per epiweek")
    args = parser.parse_args()


    input1 = args.genome_matrix
    input2 = args.case_matrix
    unique_id = args.index_column
    output1 = args.output1
    output2 = args.output2
    output3 = args.output3
    baseline = args.baseline


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_samplingbias/subsampling/run4_mundo/'
    # input1 = path + 'matrix_genomes_epiweeks.tsv'
    # input2 = path + 'matrix_cases_epiweeks.tsv'
    # output1 = path + 'weekly_sampling_proportions.tsv'
    # output2 = path + 'weekly_sampling_bias.tsv'
    # output3 = path + 'matrix_genomes_epiweeks_corrected.tsv'
    # unique_id = 'iso'
    # baseline = 0.01


    # input genome and case counts per epiweek
    separator = '\t'
    dfG = pd.read_csv(input1, encoding='utf-8', sep='\t', dtype=str)
    dfC = pd.read_csv(input2, encoding='utf-8', sep='\t', dtype=str)


    # get total genomes and cases
    date_intersection = []
    for column in dfG.columns.to_list():
        if column[-1].isdecimal():
            if column in dfC.columns.to_list():
                date_intersection.append(column)
    # print(date_intersection)

    def get_sum(df):
        df = df[date_intersection]
        df = df.astype(int)
        return df.values.sum()


    # calculate average sampling proportion
    global_samp_prop = get_sum(dfG)/get_sum(dfC) # genomes divided by cases

    # consider user defined baseline sampling proportion
    if baseline != None:
        global_samp_prop = baseline

    print('\n### Target sampling proportion:\n\n - ' + str(global_samp_prop) + '\n')

    # set new index
    dfG.set_index(unique_id, inplace=True)
    dfC.set_index(unique_id, inplace=True)


    nonDateCols = [column for column in dfG.columns.to_list() if not column[-1].isdecimal()]
    # datecols = [column for column in dfG.columns.to_list() if column[-1].isdecimal()]

    # create empty dataframes
    dfP = dfG.filter(nonDateCols, axis=1) # sampling proportion dataframe
    dfB = dfG.filter(nonDateCols, axis=1) # sampling bais dataframe
    dfW = dfG.filter(nonDateCols, axis=1) # corrected genome count dataframe

    # print(dfP)
    # print(dfB)

    # get sampling proportions and biases
    no_casedata = []
    for idx, row in dfG.iterrows():
        # for epiweek in time_cols:
        total_genomes = 0
        total_cases = 0
        for epiweek in date_intersection:
            # print(idx)
            genome_count = int(dfG.loc[idx, epiweek])
            try:
                case_count = int(dfC.loc[idx, epiweek])
            except:
                case_count = 0
                # print(idx)
                if idx not in no_casedata:
                    no_casedata.append(idx)

            samp_prop = ''
            bias = ''
            corrected_count = ''

            if int(case_count) > 0 and int(genome_count) > 0:
                if int(genome_count) > int(case_count):
                    case_count = genome_count

                samp_prop = int(genome_count)/int(case_count)
                # samp_prop = np.absolute(np.log(int(genome_count)/int(case_count)))
                bias = float(samp_prop - global_samp_prop)
                corrected_count = int(np.ceil(case_count * global_samp_prop))
                # print(genome_count, case_count, samp_prop)
                # print(idx, bias)
            elif int(case_count) > 0 and int(genome_count) == 0:
                samp_prop = 0
                bias = '-'
                corrected_count = int(np.ceil(case_count * global_samp_prop))
            else:
                samp_prop = 'X'
                bias = 'X'
                corrected_count = 0
                # print(genome_count, case_count, samp_prop)

            dfP.loc[idx, epiweek] = samp_prop # add observed sampling proportion
            dfB.loc[idx, epiweek] = bias # add calculated sampling bias

            dfW.loc[idx, epiweek] = corrected_count # add corrected genome count
            # print(corrected_count)
            dfW[epiweek] = pd.to_numeric(dfW[epiweek], downcast='integer', errors='ignore')

            # get total counts
            total_genomes += genome_count
            total_cases += case_count
        if total_cases > 0:
            dfP.loc[idx, 'cumulative_proportion'] = total_genomes / total_cases
            # print(total_genomes / total_cases)
        else:
            dfP.loc[idx, 'cumulative_proportion'] = 'NA'

    # print(dfP)
    # print(dfB)
    # print(dfW)


    # output processed dataframes
    dfP.to_csv(output1, sep='\t', index=True)
    dfB.to_csv(output2, sep='\t', index=True)
    dfW.to_csv(output3, sep='\t', index=True)

    # report
    if len(no_casedata) > 0:
        print('\n### Not case data found for:\n')
        [print(' - ' + loc) for loc in no_casedata]

    # dfP = dfP.replace('X', np.NaN).replace(0, np.NaN)
    # dfT = dfP[date_intersection].stack()
    # mean_prop = dfT.values.mean()
    # std_prop = dfT.values.std(ddof=0)
    #
    # print(mean_prop, std_prop)
    #
    # dfZ = dfG.filter(nonDateCols, axis=1) # corrected genome count dataframe
    #
    # for epiweek in sorted(date_intersection):
    #     for idx, row in dfP.iterrows():
    #         prop = dfP.loc[idx, epiweek]
    #         if prop > 0:
    #             print(prop)
    #             zscore = (prop - mean_prop) / std_prop
    #             dfZ.loc[idx, epiweek] = zscore # add observed sampling proportion
    #             # print(prop, zscore, mean_prop, std_prop)
    #             print(zscore)
    #
    #
    # dfZ.to_csv(output2, sep='\t', index=True)
