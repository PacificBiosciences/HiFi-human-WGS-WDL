version 1.0

import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "sample_analysis/sample_analysis.wdl" as SampleAnalysis
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis
import "tertiary_analysis/tertiary_analysis.wdl" as TertiaryAnalysis

workflow humanwgs {
	input {
		Cohort cohort

		ReferenceData reference
		SlivarData? slivar_data

		String deepvariant_version = "1.5.0"
		DeepVariantModel? deepvariant_model

		Int? pbsv_call_mem_gb
		Int? glnexus_mem_gb

		Boolean run_tertiary_analysis = false

		# Backend configuration
		String backend
		String? zones
		String? aws_spot_queue_arn
		String? aws_on_demand_queue_arn
		String? container_registry

		Boolean preemptible
	}

	call BackendConfiguration.backend_configuration {
		input:
			backend = backend,
			zones = zones,
			aws_spot_queue_arn = aws_spot_queue_arn,
			aws_on_demand_queue_arn = aws_on_demand_queue_arn,
			container_registry = container_registry
	}

	RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

	scatter (sample in cohort.samples) {
		call SampleAnalysis.sample_analysis {
			input:
				sample = sample,
				reference = reference,
				deepvariant_version = deepvariant_version,
				deepvariant_model = deepvariant_model,
				default_runtime_attributes = default_runtime_attributes
		}
	}

	if (length(cohort.samples) > 1) {

		scatter (sample in cohort.samples) {
			String sample_id = sample.sample_id
		}
		
		call CohortAnalysis.cohort_analysis {
			input:
				cohort_id = cohort.cohort_id,
				sample_ids = sample_id,
				aligned_bams = flatten(sample_analysis.aligned_bams),
				svsigs = flatten(sample_analysis.svsigs),
				gvcfs = sample_analysis.small_variant_gvcf,
				reference = reference,
				pbsv_call_mem_gb = pbsv_call_mem_gb,
				glnexus_mem_gb = glnexus_mem_gb,
				default_runtime_attributes = default_runtime_attributes
		}
	}

	if (run_tertiary_analysis && defined(slivar_data) && defined(reference.gnomad_af) && defined(reference.hprc_af) && defined(reference.gff) && defined(reference.population_vcfs)) {
		IndexData slivar_small_variant_input_vcf = select_first([
			cohort_analysis.phased_joint_small_variant_vcf,
			sample_analysis.phased_small_variant_vcf[0]
		])
		IndexData slivar_sv_input_vcf = select_first([
			cohort_analysis.phased_joint_sv_vcf,
			sample_analysis.phased_sv_vcf[0]
		])

		call TertiaryAnalysis.tertiary_analysis {
			input:
				cohort = cohort,
				small_variant_vcf = slivar_small_variant_input_vcf,
				sv_vcf = slivar_sv_input_vcf,
				reference = reference,
				slivar_data = select_first([slivar_data]),
				default_runtime_attributes = default_runtime_attributes
		}
	}

	output {
		# sample_analysis output

		# per movie stats, alignments
		Array[Array[File]] bam_stats = sample_analysis.bam_stats
		Array[Array[File]] read_length_summary = sample_analysis.read_length_summary
		Array[Array[File]] read_quality_summary = sample_analysis.read_quality_summary

		# per sample small variant calls
		Array[IndexData] small_variant_gvcfs = sample_analysis.small_variant_gvcf
		Array[File] small_variant_vcf_stats = sample_analysis.small_variant_vcf_stats
		Array[File] small_variant_roh_out = sample_analysis.small_variant_roh_out
		Array[File] small_variant_roh_bed = sample_analysis.small_variant_roh_bed

		# per sample final phased variant calls and haplotagged alignments
		Array[IndexData] sample_phased_small_variant_vcfs = sample_analysis.phased_small_variant_vcf
		Array[IndexData] sample_phased_sv_vcfs = sample_analysis.phased_sv_vcf
		Array[File] sample_hiphase_stats = sample_analysis.hiphase_stats
		Array[File] sample_hiphase_blocks = sample_analysis.hiphase_blocks
		Array[File] sample_hiphase_haplotags = sample_analysis.hiphase_haplotags
		Array[IndexData] merged_haplotagged_bam = sample_analysis.merged_haplotagged_bam
		Array[File] haplotagged_bam_mosdepth_summary = sample_analysis.haplotagged_bam_mosdepth_summary
		Array[File] haplotagged_bam_mosdepth_region_bed = sample_analysis.haplotagged_bam_mosdepth_region_bed
		
		# per sample trgt outputs
		Array[IndexData] trgt_spanning_reads = sample_analysis.trgt_spanning_reads
		Array[IndexData] trgt_repeat_vcf = sample_analysis.trgt_repeat_vcf
		Array[File] trgt_dropouts = sample_analysis.trgt_dropouts

		# per sample cpg outputs
		Array[Array[File]] cpg_pileup_beds = sample_analysis.cpg_pileup_beds
		Array[Array[File]] cpg_pileup_bigwigs = sample_analysis.cpg_pileup_bigwigs

		# per sample paraphase outputs
		Array[File] paraphase_output_jsons = sample_analysis.paraphase_output_json
		Array[IndexData] paraphase_realigned_bams = sample_analysis.paraphase_realigned_bam
		Array[Array[File]] paraphase_vcfs = sample_analysis.paraphase_vcfs

		# per sample hificnv outputs
		Array[IndexData] hificnv_vcfs = sample_analysis.hificnv_vcf
		Array[File] hificnv_copynum_bedgraphs = sample_analysis.hificnv_copynum_bedgraph
		Array[File] hificnv_depth_bws = sample_analysis.hificnv_depth_bw
		Array[File] hificnv_maf_bws = sample_analysis.hificnv_maf_bw

		# cohort_analysis output
		IndexData? cohort_sv_vcf = cohort_analysis.phased_joint_sv_vcf
		IndexData? cohort_small_variant_vcf = cohort_analysis.phased_joint_small_variant_vcf
		File? cohort_hiphase_stats = cohort_analysis.hiphase_stats
		File? cohort_hiphase_blocks = cohort_analysis.hiphase_blocks

		# tertiary_analysis output
		IndexData? filtered_small_variant_vcf = tertiary_analysis.filtered_small_variant_vcf
		IndexData? compound_het_small_variant_vcf = tertiary_analysis.compound_het_small_variant_vcf
		File? filtered_small_variant_tsv = tertiary_analysis.filtered_small_variant_tsv
		File? compound_het_small_variant_tsv = tertiary_analysis.compound_het_small_variant_tsv
		IndexData? filtered_svpack_vcf = tertiary_analysis.filtered_svpack_vcf
		File? filtered_svpack_tsv = tertiary_analysis.filtered_svpack_tsv
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		reference: {help: "Reference genome data"}
		slivar_data: {help: "Data files used for annotation with slivar (required if `run_tertiary_analysis` is set to `true`)"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		pbsv_call_mem_gb: {help: "Optional amount of RAM in GB for pbsv_call; default 64 for cohorts N<=3, 96 for cohorts N>3"}
		glnexus_mem_gb: {help: "Optional amount of RAM in GB for glnexus; default 30"}
		run_tertiary_analysis: {help: "Run the optional tertiary analysis steps"}
		backend: {help: "Backend where the workflow will be executed ['GCP', 'Azure', 'AWS', 'HPC']"}
		zones: {help: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"}
		aws_spot_queue_arn: {help: "Queue ARN for the spot batch queue; required if backend is set to 'AWS'"}
		aws_on_demand_queue_arn: {help: "Queue ARN for the on demand batch queue; required if backend is set to 'AWS'"}
		container_registry: {help: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used."}
		preemptible: {help: "Where possible, run tasks preemptibly"}
	}
}
