![alt text](https://github.com/andersonbrito/subsampler/blob/master/images/logo.png "subsampler")

# subsampler

A tool to subsample genomic data guided by epidemiological time series data.

# Citing

If you use this tool in a publication, please cite our paper:

[Early introductions and transmission of SARS-CoV-2 variant B.1.1.7 in the United States](https://www.sciencedirect.com/science/article/pii/S0092867421004347)

# Requirements
`subsampler` runs on MacOS and Linux.


* conda
* Fasta file containing all sampled genomes
* Metadata file containing at least sample names ("`strain`"), `date`, and geographic locations (`country`, `division`, etc).
* Matrix of daily case counts per geographic location (as listed in the metadata)


# Installation
```
git clone https://github.com/andersonbrito/subsampler.git
cd config
conda env create -f subsampler.yml
conda activate subsampler
```


# Pipeline overview

![alt text](https://github.com/andersonbrito/subsampler/blob/master/images/workflow.png "subsampler")
__Figure 1. Workflow Overview__ 


## Creating case count matrix

_`subsampler` can perform subsampling using epidemiological data from any geographical level (per country, per states, etc) provided daily case counts are available_

* Read daily case data file
* Convert date format to YYYY-MM-DD
* Generate matrix of case counts, locations versus days

## Creating genome matrix

* Read genomic metadata file
* Convert date format to YYYY-MM-DD
* Generate matrix of genome counts, locations versus days


## Aggregating genomic and epidemiological data per epiweek

* Combine genomic and case counts per epidemiological week
* Drop data from time periods outside the boundaries defined by `start_date` and `end_date`.

## Correcting genomic sampling bias

* Read matrices of epiweek genomic and case counts
* Generate matrix reporting the observed sampling proportions per epiweek
* Generate matrix reporting the sampling bias (under- and oversampling) given the baseline
* Generate matrix with the corrected genome count per week, given the pre-defined baseline sampling proportion


## Perform subsampling

* Read sequence, metadata and corrected genomic count matrix
* Read lists of genomes to be kept or remove in all instances (if provided)
* Read batch removal file, to exclude genomes from certain metadata categories
* Perform subsampling guided by case counts per epiweek
* Generate subsampled sequence, and metadata file
* Generate report with number of sampled genomes per location


# Execution

To run this pipeline, users need to provide a TSV file of daily case counts similar to the format below:

**US case counts**

|code|state         |2021-01-01|2021-01-02|2021-01-03|2021-01-04|2021-01-05|...|
|----|--------------|----------|----------|----------|----------|----------|---|
|AK  |Alaska        |5         |802       |297       |264       |200       |...|
|AL  |Alabama       |4521      |3711      |2476      |2161      |5498      |...|
|AR  |Arkansas      |4304      |2000      |2033      |1306      |4107      |...|
|AS  |American Samoa|0         |0         |0         |0         |0         |...|
|AZ  |Arizona       |10060     |8883      |17234     |5158      |5932      |...|
|CA  |California    |39425     |50222     |37016     |38256     |38962     |...|
|CO  |Colorado      |3064      |2011      |2078      |2185      |3458      |...|
|CT  |Connecticut   |0         |4412      |0         |4516      |2332      |...|
|DC  |District of Columbia|269       |257       |255       |140       |262       |...|
|... |...           |...       |...       |...       |...       |...       |...|


**Global case counts**

|code|country       |2021-01-01|2021-01-02|2021-01-03|2021-01-04|2021-01-05|...|
|---|--------------|----------|----------|----------|----------|----------|---|
|ABW|Aruba         |20        |23        |32        |42        |110       |...|
|AFG|Afghanistan   |0         |0         |0         |1485      |94        |...|
|AGO|Angola        |15        |40        |34        |42        |72        |...|
|AIA|Anguilla      |0         |0         |2         |0         |0         |...|
|ALB|Albania       |0         |675       |447       |185       |660       |...|
|AND|Andorra       |68        |49        |26        |57        |59        |...|
|ARE|United Arab Emirates|1856      |1963      |1590      |1501      |1967      |...|
|ARG|Argentina     |4080      |5240      |5884      |8222      |13790     |...|
|ARM|Armenia       |329       |60        |229       |193       |324       |...|
|...|...           |...       |...       |...       |...       |...       |...|


Using one of the commands below, users can download reformatted daily case count files automatically from [CSSE at Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19):

**Global data**
```
python scripts/get_daily_matrix_global.py --download yes
```

**US data**
```
python scripts/get_daily_matrix_usa.py --download yes
```

Users can provide their own daily case count file, as long as it matches the format above (tab-separated, with daily counts, and a column with unique identifiers). If one of the commands above is used, the reformatted matrix of case counts need to be placed inside `/data`.

Now, edit the Snakefile to fix the followin lines:

* [index_column](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L10) = "code" (this should match the index column with unique identifiers)

* [end_date](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L17) = "YYYY-MM-DD" (select an end date according to the case data file)

* [extra_columns](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L32) = "second column with identifier" (one can select another column to be display alongside the `index_column`)


## Obtaining the percentage of sequenced cases per week

The `subsampler` pipeline allows users to calculate the percentage of sequenced cases in countries and US states. It aggregates both genome counts and case counts per week per location (country or state), and proceed with the division genomes/cases to get a time series of proportion of sequenced genomes, information useful for monitoring how genomic surveillance is going in different regions.


To that end, the user needs to provide a metadata matrix, similar to the one used by [nextstrain](http://nextstrain.org), which can be downloaded from [GISAID](http://gisaid.org), under `Downloads > Genomic Epidemiology`. Rename such file as `metadata_nextstrain.tsv`, place it inside `/data`, and the pipeline only half-way through, by using the command:

```
conda activate subsampler
snakemake correct_bias
```

After a few minutes, among the files in `/outputs`, users will find three matrices, one of them showing the weekly proportion of sequenced cases:

```
matrix_cases_epiweeks.tsv
matrix_genomes_epiweeks.tsv
weekly_sampling_proportions.tsv
```

## Obtaining a list of genomes, sampled based on time series of COVID-19 cases

To run the full pipeline, and obtain the list of genomes sampled based on COVID-19 case counts, the last step of the pipeline need to be executed:


```
snakemake subsample
```

