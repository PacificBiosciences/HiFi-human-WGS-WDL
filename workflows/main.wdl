version 1.0

import "common/structs.wdl"
import "sample_analysis/sample_analysis.wdl" as SampleAnalysis
import "de_novo_assembly_sample/de_novo_assembly_sample.wdl" as DeNovoAssemblySample
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis
import "de_novo_assembly_trio/de_novo_assembly_trio.wdl" as DeNovoAssemblyTrio
import "slivar/slivar.wdl" as Slivar

workflow humanwgs {
	input {
		Cohort cohort

		ReferenceData reference
		SlivarData slivar_data

		String deepvariant_version
		DeepVariantModel? deepvariant_model

		String container_registry
	}

	scatter (sample in cohort.samples) {
		call SampleAnalysis.sample_analysis {
			input:
				sample = sample,
				reference = reference,
				deepvariant_version = deepvariant_version,
				deepvariant_model = deepvariant_model,
				container_registry = container_registry
		}

		if (sample.run_de_novo_assembly) {
			call DeNovoAssemblySample.de_novo_assembly_sample {
				input:
					sample = sample,
					reference = reference,
					container_registry = container_registry
			}
		}
	}

	if (length(cohort.samples) > 1) {
		call CohortAnalysis.cohort_analysis {
			input:
				cohort_id = cohort.cohort_id,
				aligned_bams = flatten(sample_analysis.aligned_bams),
				svsigs = flatten(sample_analysis.svsigs),
				gvcfs = sample_analysis.small_variant_gvcf,
				reference = reference,
				container_registry = container_registry
		}

		if (cohort.run_de_novo_assembly_trio) {
			call DeNovoAssemblyTrio.de_novo_assembly_trio {
				input:
					cohort = cohort,
					reference = reference,
					container_registry = container_registry
			}
		}
	}

	IndexData slivar_small_variant_input_vcf = select_first([
		cohort_analysis.phased_joint_called_vcf,
		sample_analysis.phased_small_variant_vcf[0]
	])
	IndexData slivar_sv_input_vcf = select_first([
		cohort_analysis.sv_vcf,
		sample_analysis.sv_vcf[0]
	])

	call Slivar.slivar {
		input:
			cohort = cohort,
			small_variant_vcf = slivar_small_variant_input_vcf,
			sv_vcf = slivar_sv_input_vcf,
			reference = reference,
			slivar_data = slivar_data,
			container_registry = container_registry
	}

	output {
		# sample_analysis.smrtcells_analysis output
		Array[Array[File]] bam_stats = sample_analysis.bam_stats
		Array[Array[File]] read_length_summary = sample_analysis.read_length_summary
		Array[Array[File]] read_quality_summary = sample_analysis.read_quality_summary
		Array[Array[IndexData]] aligned_bams = sample_analysis.aligned_bams
		Array[Array[File]] aligned_bam_mosdepth_summary = sample_analysis.aligned_bam_mosdepth_summary
		Array[Array[File]] aligned_bam_mosdepth_region_bed = sample_analysis.aligned_bam_mosdepth_region_bed
		Array[IndexData] small_variant_vcfs = sample_analysis.small_variant_vcf
		Array[IndexData] small_variant_gvcfs = sample_analysis.small_variant_gvcf
		Array[File] small_variant_vcf_stats = sample_analysis.small_variant_vcf_stats
		Array[File] small_variant_roh_bed = sample_analysis.small_variant_roh_bed

		# sample_analysis output
		Array[IndexData] sample_sv_vcfs = sample_analysis.sv_vcf
		Array[IndexData] sample_phased_small_variant_vcfs = sample_analysis.phased_small_variant_vcf
		Array[File] sample_whatshap_stats_gtfs = sample_analysis.whatshap_stats_gtf
		Array[File] sample_whatshap_stats_tsvs = sample_analysis.whatshap_stats_tsv
		Array[File] sample_whatshap_stats_blocklists = sample_analysis.whatshap_stats_blocklist
		Array[IndexData] merged_haplotagged_bam = sample_analysis.merged_haplotagged_bam
		Array[File] haplotagged_bam_mosdepth_summary = sample_analysis.haplotagged_bam_mosdepth_summary
		Array[File] haplotagged_bam_mosdepth_region_bed = sample_analysis.haplotagged_bam_mosdepth_region_bed
		Array[IndexData] trgt_spanning_reads = sample_analysis.trgt_spanning_reads
		Array[IndexData] trgt_repeat_vcf = sample_analysis.trgt_repeat_vcf
		Array[File] trgt_dropouts = sample_analysis.trgt_dropouts
		Array[Array[File]] cpg_pileups = sample_analysis.cpg_pileups

		# de_novo_assembly_sample output
		Array[Array[File]?] assembly_noseq_gfas = de_novo_assembly_sample.assembly_noseq_gfas
		Array[Array[File]?] assembly_lowQ_beds = de_novo_assembly_sample.assembly_lowQ_beds
		Array[Array[File]?] zipped_assembly_fastas = de_novo_assembly_sample.zipped_assembly_fastas
		Array[Array[File]?] assembly_stats = de_novo_assembly_sample.assembly_stats
		Array[IndexData?] asm_bam = de_novo_assembly_sample.asm_bam
		Array[IndexData?] htsbox_vcf = de_novo_assembly_sample.htsbox_vcf
		Array[File?] htsbox_vcf_stats = de_novo_assembly_sample.htsbox_vcf_stats

		# cohort_analysis output
		IndexData? cohort_sv_vcf = cohort_analysis.sv_vcf
		IndexData? cohort_phased_joint_called_vcf = cohort_analysis.phased_joint_called_vcf
		File? cohort_whatshap_stats_gtfs = cohort_analysis.whatshap_stats_gtf
		File? cohort_whatshap_stats_tsvs = cohort_analysis.whatshap_stats_tsv
		File? cohort_whatshap_stats_blocklists = cohort_analysis.whatshap_stats_blocklist

		# de_novo_assembly_trio output
		Map[String, String]? haplotype_key = de_novo_assembly_trio.haplotype_key
		Array[File]? trio_assembly_noseq_gfas = de_novo_assembly_trio.assembly_noseq_gfas
		Array[File]? trio_assembly_lowQ_beds = de_novo_assembly_trio.assembly_lowQ_beds
		Array[File]? trio_zipped_assembly_fastas = de_novo_assembly_trio.zipped_assembly_fastas
		Array[File]? trio_assembly_stats = de_novo_assembly_trio.assembly_stats
		IndexData? trio_asm_bam = de_novo_assembly_trio.asm_bam
		Array[File]? yak_trioeval = de_novo_assembly_trio.trioeval

		# slivar output
		IndexData filtered_small_variant_vcf = slivar.filtered_small_variant_vcf
		IndexData compound_het_small_variant_vcf = slivar.compound_het_small_variant_vcf
		File filtered_small_variant_tsv = slivar.filtered_small_variant_tsv
		File compound_het_small_variant_tsv = slivar.compound_het_small_variant_tsv
		IndexData filtered_svpack_vcf = slivar.filtered_svpack_vcf
		File filtered_svpack_tsv = slivar.filtered_svpack_tsv
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		reference: {help: "Reference genome data"}
		slivar_data: {help: "Data files used for annotation with slivar"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}
