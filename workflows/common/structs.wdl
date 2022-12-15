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
	Array[File] movie_bams

	String sex
	Boolean affected

	String? father_id
	String? mother_id

	Boolean run_de_novo_assembly
}

struct Cohort {
	String cohort_id
	Array[Sample] samples

	Array[String] phenotypes

	Boolean run_de_novo_assembly_trio
}

struct FamilySampleIndices {
	Array[Int] child_indices
	Int father_index
	Int mother_index
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

struct DeepVariantModel {
	IndexData model
	File metadata
}
