version 1.0

import "wdl-common/wdl/structs.wdl"

struct ReferenceData {
	String name
	IndexData fasta

	File pbsv_splits

	File tandem_repeat_bed
	File trgt_tandem_repeat_bed

	IndexData hificnv_exclude_bed
	File hificnv_expected_bed_male
	File hificnv_expected_bed_female

	File? gnomad_af
	File? hprc_af
	File? gff

	Array[IndexData]? population_vcfs
}

struct Sample {
	String sample_id
	Array[File] movie_bams

	String? sex
	Boolean affected

	String? father_id
	String? mother_id
}

struct Cohort {
	String cohort_id
	Array[Sample] samples

	Array[String] phenotypes
}

struct SlivarData {
	File slivar_js
	File hpo_terms
	File hpo_dag
	File hpo_annotations
	File ensembl_to_hgnc
	File lof_lookup
	File clinvar_lookup
}
