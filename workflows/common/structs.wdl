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

	IndexData eee_vcf
	IndexData gnomad_sv_vcf
	IndexData hprc_pbsv_vcf
	IndexData decode_vcf
}

struct Sample {
	String sample_id
	Array[IndexData] movie_bams

	String sex
	Boolean affected

	String? father_id
	String? mother_id
}

struct Cohort {
	String cohort_id
	Array[Sample] samples

	# TODO non-optional?
	Array[String]? phenotypes
}

struct HpoData {
	File terms
	File dag
	File annotations
}
