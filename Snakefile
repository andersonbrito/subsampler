rule arguments:
	params:
		sequences = "data/gisaid_hcov-19.fasta",
		metadata = "data/metadata_nextstrain.tsv",
		case_data = "data/time_series_covid19_confirmed_US_reformatted.tsv",
		include = "config/keep.txt",
		exclude = "config/remove.txt",
		drop = "config/batch_removal.tsv",
		index_column = "state_code",
		date_column = "date",
		baseline = "0.002",
		refgenome_size = "29420",
		max_missing = "5",
		seed_num = "2007",
		start_date = "2019-12-15",
		end_date = "2020-08-29"


arguments = rules.arguments.params


rule genome_matrix:
	message:
		"""
		Generate matrix of genome counts per day, for each element in column="{arguments.index_column}"
		"""
	input:
		metadata = arguments.metadata
	params:
		index = arguments.index_column,
		extra_columns = "region country",
		date = arguments.date_column
	output:
		matrix = "outputs/genome_matrix_days.tsv"
	shell:
		"""
		python3 scripts/get_genome_matrix.py \
			--metadata {input.metadata} \
			--index-column {params.index} \
			--extra-columns {params.extra_columns} \
			--date-column {params.date} \
			--output {output.matrix}
		"""


rule epiweek_conversion:
	message:
		"""
		Generate matrix of genome and case counts per epiweek
		"""
	input:
		genome_matrix = "outputs/genome_matrix_days.tsv",
		case_matrix = arguments.case_data
	output:
		output1 = "outputs/matrix_genomes_epiweeks.tsv",
		output2 = "outputs/matrix_cases_epiweeks.tsv"
	params:
		start_date = "2020-02-22"
	shell:
		"""
		python3 scripts/convertDays2Epiweek.py \
			--input {input.genome_matrix} \
			--output {output.output1}
		python3 scripts/convertDays2Epiweek.py \
			--input {input.case_matrix} \
			--start-date {params.start_date} \
			--output {output.output2}
		"""


rule correct_bias:
	message:
		"""
		Correct under- and oversampling genome counts based on epidemiological data
		"""
	input:
		genome_matrix = "outputs/matrix_genomes_epiweeks.tsv",
		case_matrix = "outputs/matrix_cases_epiweeks.tsv"
	params:
		index = arguments.index_column,
		baseline = arguments.baseline
	output:
		output1 = "outputs/weekly_sampling_proportions.tsv",
		output2 = "outputs/weekly_sampling_bias.tsv",
		output3 = "outputs/matrix_genomes_epiweeks_corrected.tsv"
	shell:
		"""
		python3 scripts/correct_bias.py \
			--genome-matrix {input.genome_matrix} \
			--case-matrix {input.case_matrix} \
			--index-column {params.index} \
			--baseline {params.baseline} \
			--output1 {output.output1} \
			--output2 {output.output2} \
			--output3 {output.output3}
		"""


rule subsample:
	message:
		"""
		Sample genomes and metadata according to the corrected genome matrix
		"""
	input:
		sequences = arguments.sequences,
		metadata = arguments.metadata,
		corrected_matrix = "outputs/matrix_genomes_epiweeks_corrected.tsv",
		keep = arguments.include,
		remove = arguments.exclude
	params:
		size = arguments.refgenome_size,
		missing = arguments.max_missing,
		seed = arguments.seed_num,
		index = arguments.index_column,
		date = arguments.date_column,
		filter = arguments.date_column, # change it if another date column needs to be used
		start = arguments.start_date,
		end = arguments.end_date
	output:
		output1 = "outputs/sequences.fasta",
		output2 = "outputs/metadata.tsv",
		output3 = "outputs/sampling_stats.txt"
	shell:
		"""
		python3 scripts/subsampler_timeseries.py \
			--sequences {input.sequences} \
			--metadata {input.metadata} \
			--genome-matrix {input.corrected_matrix} \
			--max-missing {params.missing} \
			--refgenome-size {params.size} \
			--keep {input.keep} \
			--remove {input.remove} \
			--seed {params.seed} \
			--index-column {params.index} \
			--date-column {params.date} \
			--filter-column {params.filter} \
			--start-date {params.start} \
			--end-date {params.end} \
			--sampled-sequences {output.output1} \
			--sampled-metadata {output.output2} \
			--report {output.output3}
		echo '# Sampling proportion: {arguments.baseline}' | cat - {output.output3} > temp && mv temp {output.output3}
		"""


rule clean:
	message: "Removing directories: {params}"
	params:
		"outputs"
	shell:
		"""
		rm -rfv {params}
		"""
