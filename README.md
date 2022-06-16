![alt text](https://github.com/andersonbrito/subsampler/blob/master/images/logo.png "subsampler")

# subsampler

A pipeline for subsampling genomic data based on epidemiological time series data.

# Citing

If you use this tool in a publication, please cite our paper:
> Alpert, T., Brito, A. F., Lasek-Nesselquist, E., Rothman, J., Valesano, A. L., MacKay, M. J., ... & Grubaugh, N. D. (2021). [Early introductions and transmission of SARS-CoV-2 variant B.1.1.7 in the United States](https://www.sciencedirect.com/science/article/pii/S0092867421004347). Cell, 184(10), 2595-2604.


# Requirements
`subsampler` runs on MacOS and Linux. To run all steps until the actual subsampling, besides having `conda` and the `subsampler` environment installed (see next section), you need to provide the following files:

* Metadata file containing columns with at least: sample names (`strain`, or `accession number`), `date`, and geographic columns (`country`, `division`, etc).
* Matrix of daily case counts per geographic location (matching the geographic level of interested, included in the metadata)


---
## Note

1. If you only need to run this pipeline to calculate the proportion of sequenced cases per geographic location, per unit of time, you just need to run the pipeline up to the `correct_bias` step (`snakemake correct_bias`). It will produce a matrix with the proportions of sequenced cases.

2. For this pipeline to run up to its last step (`snakemake subsample`), the metadata file itself must contain the minimum set of columns (described above), correctly named in the Snakefile ([here](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L9)).

3. The complete run of `subsampler` generates, among other files, a TXT file containing a list of accession numbers (`gisaid_epi_isl`, for example) or genome names (`strain`, for example), provided the corresponding columns are found in the metadata file used as input. Such list of genome names or (specially) accession numbers can be used to download an actual sequence file from a genomic database (GISAID, for example), so that further analyses can be performed.


---

# Installation
```
git clone https://github.com/andersonbrito/subsampler.git
cd config
conda env create -f subsampler.yaml
conda activate subsampler
```

To update the conda environment, run:

```
cd config
conda env update -f subsampler.yaml
```

Alternatively, `mamba` can also be used to install the `subsampler` conda environment.

# Pipeline overview

![alt text](https://github.com/andersonbrito/subsampler/blob/master/images/workflow.png "subsampler")
__Figure 1. Workflow Overview__ 


## Creating case count matrix

_`subsampler` can perform subsampling using epidemiological data from any geographical level (per country, per states, etc) provided daily case counts are available_. See more details in the 'Execution' section. To prepare

* Download and provide a daily case data file
* Generate matrix of case counts, per location (Y axis), per day (X axis)

## Creating genome matrix

The pipeline will perform these actions:
* Read genomic metadata file
* Convert date format to YYYY-MM-DD
* Generate matrix of genome counts, per location (Y axis), per day (X axis)


## Aggregating genomic and epidemiological data per epiweek

The pipeline will perform these actions:
* Combine genomic and case counts per unit of time (week, month, or year)
* Drop data from time periods outside the boundaries defined by `start_date` and `end_date`.

## Correcting genomic sampling bias

The pipeline will perform these actions:
* Read matrices with case and genome counts
* Generate matrix reporting the observed sampling proportions per unit of time
* Generate matrix reporting the sampling bias (under- and oversampling) given the baseline defined by the user
* Generate matrix with the corrected genome count per unit of time, given the pre-defined baseline sampling proportion


## Perform subsampling

The pipeline will perform these actions:
* Read metadata and corrected genomic count matrix
* Read lists of genomes to be kept or remove in all instances (if provided)
* Read the filter file, to include or exclude genomes from certain metadata categories
* Perform subsampling guided by case counts per unit of time
* Generate a list subsampled sequences, and a corresponding metadata file
* Generate a report with number of sampled genomes per location

# Execution

To run this pipeline, users need to provide a TSV file of daily case counts similar to the format below:

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


**Country-level case counts**

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



Using one of the commands below, users can download reformatted daily case count files automatically from [CSSE at Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19):

**Global data**
```
python scripts/get_daily_matrix_global.py --download yes
```

**US data**
```
python scripts/get_daily_matrix_usa.py --download yes
```

Users can provide their own daily case count file, as long as it matches the format above (tab-separated, with daily counts, and a column with unique identifiers). If one of the commands above is used, the reformatted matrix of case counts would need to be placed inside `/data`, and should be named [here](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L5).

Now, edit the Snakefile to fix the following lines:

* [start_date](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L16) = "YYYY-MM-DD" (select the start date according to your needs)

* [end_date](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L17) = "YYYY-MM-DD" (select the end date according to your needs)

* [extra_columns](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L32) = second column with identifier, such as region, continent (a column found in the original metadata file, which you want to see displayed alongside the `geo_column` in the final outputs)


## Obtaining the percentage of sequenced cases per week

The `subsampler` pipeline allows users to calculate the percentage of sequenced cases per location. It aggregates both genome counts and case counts per unit of time, per location (country or state), and proceed with calculations (genomes/cases) to get a time series of proportions of sequenced cases, information useful for monitoring how genomic surveillance is going in different locations.


To that end, the user needs to provide a metadata matrix, similar to the one used by [nextstrain](http://nextstrain.org), which can be downloaded from [GISAID](http://gisaid.org), under `Downloads > Genomic Epidemiology`. Add the name of that file [here](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L4), place it inside `/data`, and run the pipeline up to the `correct_bias` step using the command below:

```
snakemake correct_bias
```

After a few minutes, among the files in `/outputs`, users will find three matrices, one of them showing the weekly proportion of sequenced cases:

```
matrix_cases_unit.tsv
matrix_genomes_unit.tsv
weekly_sampling_proportions.tsv
```

## Obtaining a list of genomes, sampled based on time series of COVID-19 cases

To obtain a list of genomes sampled based on case counts, the last step of the pipeline need to be executed:


```
snakemake subsample
```

### Downloading genome sequences

One of the outputs of this pipeline is `selected_sequences.txt`. If `id_column` is set [here](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L9) as `gisaid_epi_isl`, a list of accession numbers of subsampled genomes will be generated. Using that list, proceed as follows:

1. Go to [gisaid.org](https://www.gisaid.org/); and visit the 'Search' page of 'EpiCov'.
2. Click on 'Select'; paste the list of accession numbers in the search box; click on 'OK'; and choose the format 'Sequences (FASTA)'.

![alt text](https://github.com/andersonbrito/subsampler/blob/master/images/gisaid_download.png?raw=true)



## Subsampler may not be what you need if...

If your questions are not directly related to phylogeography, the `subsampler` approach (to obtain genomes sampled based on case counts) may not be what you need. Since the sampling is weighted by case counts, `subsampler` is more likely to sample genomes from heavily impacted countries (those with more reported cases), and the lower the [baseline](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L12) (the percentage of sequenced cases) being used, the less likely would it be for countries facing small scale outbreaks to be represented (for example, the least populated countries), which end up being overshadowed by larger countries, which may report hundreds of thousands of cases per week.

If you are not trying to infer ancestral states in a phylogeographic perspective, but is more interested in phylodynamic questions, subsampling based on the timing of the events (waves, introductions, seasons, etc) is a better approach.

In this repository you can find `genome_selector.py`, a [python script](https://github.com/andersonbrito/subsampler/blob/master/scripts/genome_selector.py) designed to sample genomes following specific variables (date of colection, country, viral lineage, or any metadata column). For example, the table below illustrates a sampling scheme to obtain around 650 genomes of viruses belonging to lineage `B.1.1.7` (alpha variant), circulating in the US and the UK, between 2020-12-01 and 2021-06-30, with other US and European samples as contextual genomes. In this hypothetical example, note that contextual genomes are selected from two time periods, and in different proportions: 50 genomes up to late November 2020, and 100 genomes from December 2020 onwards (from any lineage).

Also, the scheme is set up to ignore genomes from California and Scotland: genomes from those locations will not be included in any instance, since they are filtered out prior to the genome selection step. To reproduce the scheme below, `genome_selector.py` will use a `--metadata` file listing all genomes from the locations and lineages represented below.

`genome_selector.py` is not part of `subsampler`, and should be executed separately:

```
genome_selector.py [-h] --metadata METADATA [--keep KEEP] [--remove REMOVE] --scheme SCHEME [--report REPORT]
```
... where `--scheme` is a TSV file with this format:

|purpose|filter       |value     |filter2|value2    |sample_size|start     |end       |
|-------|-------------|----------|-------|----------|-----------|----------|----------|
|focus  |pango_lineage|B.1.1.7   |country|USA       |200        |2020-12-01|2021-06-30|
|focus  |pango_lineage|B.1.1.7   |country|United Kingdom|200        |2020-12-01|2021-06-30|
|context|region       |Europe    |       |          |100        |2020-12-01|2021-06-30|
|context|country      |USA       |       |          |50        |          |2020-11-30|
|context|country      |USA       |       |          |100        |2020-12-01|          |
|ignore |division     |California|       |          |           |          |          |
|ignore |division     |Scotland  |       |          |           |          |          |


Among the outputs of `genome_selector.py` users will find text files containing a list of around 650 genome names (e.g. USA/CT-CDC-LC0062417/2021) and another with the corresponding accession numbers (e.g. EPI_ISL_2399048). The last one can be used to filter and download genomes directly from [gisaid.org](https://www.gisaid.org/), as explained [above](https://github.com/andersonbrito/subsampler#downloading-genome-sequences).



## Latest major updates

2022-06-12:
* Fasta file with the actual sequences are no longer required as input file. Now, by default, the pipeline will not inspect the level of completeness of the genomes, but will focus on subsampling based on metadata rows only. However, asessement of sequence quality is still supported.
* A `filter_file` is now an input of this pipeline (see 'config/filters.tsv', and this [line](https://github.com/andersonbrito/subsampler/blob/master/Snakefile#L8). With this file, users can determine specific data categories to be included or excluded. This feature is useful, for example, for subsampling 'variant-specific' data (e.g. include → pango_lineage → B.1.1.7), among other uses.
* Uses can now specify the time unit of the time series (week, month or year), which should be set according to the most adequate time period that match the evolutionary time scale of the viruses under study (for SARS-CoV-2, 'week' is an adequate option, but for Dengue virus, 'month' is the best option).
