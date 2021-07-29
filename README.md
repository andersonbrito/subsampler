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

# Execution

```
conda activate subsampler
snakemake subsample
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



<!-- 
## Outputs:

#### <OUTPUT_PREFIX>.fastq

A merge of all files in the fastq directory specified as input.

#### <OUTPUT_PREFIX>_periscope_counts.csv

The counts of genomic, sub-genomic and normalisation values for known ORFs

#### <OUTPUT_PREFIX>_periscope_amplicons.csv

The amplicon by amplicon counts, this file is useful to see where the counts come from. Multiple amplicons may be represented more than once where they may have contributed to more than one ORF.

#### <OUTPUT_PREFIX>_periscope_novel_counts.csv

The counts of genomic, sub-genomic and normalisation values for non-canonical ORFs

#### <OUTPUT_PREFIX>.bam

minmap2 mapped reads and index with no adjustments made.

#### <OUTPUT_PREFIX>_periscope.bam

This is the original input bam file and index created by periscope with the reads specified in the fastq-dir. This file, however, has tags which represent the results of periscope:

- XS is the alignment score
- XA is the amplicon number
- XC is the assigned class (gDNA or sgDNA)
- XO is the orf assigned

These are useful for manual review in IGV or similar genome viewer. You can sort or colour reads by these tags to aid in manual review and figure creation.
 -->


<!-- 
# Citations

**Title**
Authors. Journal; doi;
 -->
