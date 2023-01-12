version 1.0

import "common/structs.wdl"
import "common/backend_configuration.wdl" as BackendConfiguration
import "sample_analysis/sample_analysis.wdl" as SampleAnalysis
import "de_novo_assembly_sample/de_novo_assembly_sample.wdl" as DeNovoAssemblySample
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis
import "de_novo_assembly_trio/de_novo_assembly_trio.wdl" as DeNovoAssemblyTrio
import "tertiary_analysis/tertiary_analysis.wdl" as TertiaryAnalysis

workflow humanwgs {
	input {
		Cohort cohort

		ReferenceData reference
		SlivarData slivar_data

		String deepvariant_version
		DeepVariantModel? deepvariant_model

		Int? assembly_threads

		# Backend configuration
		String backend
		String? zones
		String? aws_spot_queue_arn
		String? aws_on_demand_queue_arn
		Boolean preemptible
	}

	call BackendConfiguration.backend_configuration {
		input:
			backend = backend,
			zones = zones,
			aws_spot_queue_arn = aws_spot_queue_arn,
			aws_on_demand_queue_arn = aws_on_demand_queue_arn
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

		if (sample.run_de_novo_assembly) {
			call DeNovoAssemblySample.de_novo_assembly_sample {
				input:
					sample = sample,
					reference = reference,
					assembly_threads = assembly_threads,
					default_runtime_attributes = default_runtime_attributes,
					on_demand_runtime_attributes = backend_configuration.on_demand_runtime_attributes
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
				default_runtime_attributes = default_runtime_attributes
		}

		if (cohort.run_de_novo_assembly_trio) {
			call DeNovoAssemblyTrio.de_novo_assembly_trio {
				input:
					cohort = cohort,
					reference = reference,
					assembly_threads = assembly_threads,
					default_runtime_attributes = default_runtime_attributes,
					on_demand_runtime_attributes = backend_configuration.on_demand_runtime_attributes
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

	call TertiaryAnalysis.tertiary_analysis {
		input:
			cohort = cohort,
			small_variant_vcf = slivar_small_variant_input_vcf,
			sv_vcf = slivar_sv_input_vcf,
			reference = reference,
			slivar_data = slivar_data,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		# sample_analysis.smrtcells_analysis output
		Array[Array[File]] bam_stats = sample_analysis.bam_stats
		Array[Array[File]] read_length_summary = sample_analysis.read_length_summary
		Array[Array[File]] read_quality_summary = sample_analysis.read_quality_summary
		Array[Array[IndexData]] aligned_bams = sample_analysis.aligned_bams
		Array[Array[File]] aligned_bam_mosdepth_summary = sample_analysis.aligned_bam_mosdepth_summary
		Array[Array[File]] aligned_bam_mosdepth_region_bed = sample_analysis.aligned_bam_mosdepth_region_bed

		# sample_analysis output
		Array[IndexData] small_variant_vcfs = sample_analysis.small_variant_vcf
		Array[IndexData] small_variant_gvcfs = sample_analysis.small_variant_gvcf
		Array[File] small_variant_vcf_stats = sample_analysis.small_variant_vcf_stats
		Array[File] small_variant_roh_bed = sample_analysis.small_variant_roh_bed
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
		Array[Map[String, String]]? haplotype_key = de_novo_assembly_trio.haplotype_key
		Array[Array[File]]? trio_assembly_noseq_gfas = de_novo_assembly_trio.assembly_noseq_gfas
		Array[Array[File]]? trio_assembly_lowQ_beds = de_novo_assembly_trio.assembly_lowQ_beds
		Array[Array[File]]? trio_zipped_assembly_fastas = de_novo_assembly_trio.zipped_assembly_fastas
		Array[Array[File]]? trio_assembly_stats = de_novo_assembly_trio.assembly_stats
		Array[IndexData]? trio_asm_bams = de_novo_assembly_trio.asm_bams

		# tertiary_analysis output
		IndexData filtered_small_variant_vcf = tertiary_analysis.filtered_small_variant_vcf
		IndexData compound_het_small_variant_vcf = tertiary_analysis.compound_het_small_variant_vcf
		File filtered_small_variant_tsv = tertiary_analysis.filtered_small_variant_tsv
		File compound_het_small_variant_tsv = tertiary_analysis.compound_het_small_variant_tsv
		IndexData filtered_svpack_vcf = tertiary_analysis.filtered_svpack_vcf
		File filtered_svpack_tsv = tertiary_analysis.filtered_svpack_tsv
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		reference: {help: "Reference genome data"}
		slivar_data: {help: "Data files used for annotation with slivar"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		assembly_threads: {help: "Number of threads to use for de novo assembly"}
		backend: {help: "Backend where the workflow will be executed ['GCP', 'Azure', 'AWS']"}
		zones: {help: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"}
		aws_spot_queue_arn: {help: "Queue ARN for the spot batch queue; required if backend is set to 'AWS'"}
		aws_on_demand_queue_arn: {help: "Queue ARN for the on demand batch queue; required if backend is set to 'AWS'"}
		preemptible: {help: "Where possible, run tasks preemptibly"}
	}
}
