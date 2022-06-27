rule arguments:
	params:
#		sequences = "",
		metadata = "data/metadata.tsv",
		case_data = "data/time_series_covid19_global_reformatted.tsv",
		keep_file = "config/keep.txt",
		remove_file = "config/remove.txt",
		filter_file = "config/filters.tsv",
		id_column = "gisaid_epi_isl",
		geo_column = "country_exposure",
		date_column = "date",
		baseline = "0.001",
		refgenome_size = "1",
		max_missing = "99",
		seed_num = "2007",
		start_date = "2020-03-01",
		end_date = "2021-12-31",
		unit = "week"


arguments = rules.arguments.params


rule genome_matrix:
	message:
		"""
		Generate matrix of genome counts per day, for each element in column="{arguments.geo_column}"
		"""
	input:
		metadata = arguments.metadata
	params:
		index = arguments.geo_column,
		extra_columns = "country_exposure",
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


rule unit_conversion:
	message:
		"""
		Generate matrix of genome and case counts per epiweek
		"""
	input:
		genome_matrix = "outputs/genome_matrix_days.tsv",
		case_matrix = arguments.case_data
	output:
		output1 = "outputs/matrix_genomes_unit.tsv",
		output2 = "outputs/matrix_cases_unit.tsv"
	params:
		start_date = "2020-02-22",
		format = "integer",
		time_unit = arguments.unit
	shell:
		"""
		python3 scripts/aggregator.py \
			--input {input.genome_matrix} \
			--unit {arguments.unit} \
			--format {params.format} \
			--output {output.output1}

		python3 scripts/aggregator.py \
			--input {input.case_matrix} \
			--unit {arguments.unit} \
			--format {params.format} \
			--start-date {params.start_date} \
			--output {output.output2}
		"""


rule correct_bias:
	message:
		"""
		Correct under- and oversampling genome counts based on epidemiological data
		"""
	input:
		genome_matrix = "outputs/matrix_genomes_unit.tsv",
		case_matrix = "outputs/matrix_cases_unit.tsv"
	params:
		index = 'code',
		baseline = arguments.baseline
	output:
		output1 = "outputs/weekly_sampling_proportions.tsv",
		output2 = "outputs/weekly_sampling_bias.tsv",
		output3 = "outputs/matrix_genomes_unit_corrected.tsv"
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
#		sequences = arguments.sequences,
		metadata = arguments.metadata,
		corrected_matrix = "outputs/matrix_genomes_unit_corrected.tsv",
		keep = arguments.keep_file,
		remove = arguments.remove_file,
		filter_file = arguments.filter_file,
	params:
		size = arguments.refgenome_size,
		missing = arguments.max_missing,
		seed = arguments.seed_num,
		id_column = arguments.id_column,
		geo_column = arguments.geo_column,
		date = arguments.date_column,
		start = arguments.start_date,
		end = arguments.end_date,
		time_unit = arguments.unit,
		weekasdate = 'no'
	output:
		output1 = "outputs/selected_sequences.txt",
		output2 = "outputs/selected_metadata.tsv",
		output3 = "outputs/sampling_stats.txt"
	shell:
		"""
		python3 scripts/subsampler_timeseries.py \
			--metadata {input.metadata} \
			--genome-matrix {input.corrected_matrix} \
			--max-missing {params.missing} \
			--refgenome-size {params.size} \
			--keep {input.keep} \
			--remove {input.remove} \
			--filter-file {input.filter_file} \
			--seed {params.seed} \
			--index-column {params.id_column} \
			--geo-column {params.geo_column} \
			--date-column {params.date} \
			--time-unit {params.time_unit} \
			--weekasdate {params.weekasdate} \
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
