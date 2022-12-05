version 1.0

struct IndexData {
	File data
	File data_index
}

struct ReferenceData {
	String name
	IndexData fasta

	Array[String] chromosomes
	File chromosome_lengths

	File tandem_repeat_bed
	File trgt_tandem_repeat_bed

	File gnomad_af
	File hprc_af
	File gff
}

struct Sample {
	String sample_id
	Array[IndexData] movie_bams

	# TODO non-optional?
	Array[String]? phenotypes
	Boolean? affected

	# TODO keep this info in a PED?
	String? mother_id
	String? father_id
}

struct Cohort {
	String cohort_id
	Array[Sample] samples

	# TODO auto-generate this file based on the cohort information?
	File pedigree
}
