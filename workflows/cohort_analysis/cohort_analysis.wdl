version 1.0

# Run joint genotyping for a cohort. This workflow will be run if there is more than one sample in the cohort.

import "../humanwgs_structs.wdl"
import "../wdl-common/wdl/tasks/pbsv_call.wdl" as PbsvCall
import "../wdl-common/wdl/tasks/zip_index_vcf.wdl" as ZipIndexVcf
import "../wdl-common/wdl/tasks/glnexus.wdl" as Glnexus
import "../wdl-common/wdl/workflows/phase_vcf/phase_vcf.wdl" as PhaseVcf

workflow cohort_analysis {
	input {
		String cohort_id
		Int sample_count
		Array[IndexData] aligned_bams
		Array[File] svsigs
		Array[IndexData] gvcfs

		ReferenceData reference

		Int? pbsv_call_mem_gb
		Int? glnexus_mem_gb

		RuntimeAttributes default_runtime_attributes
	}

	scatter (gvcf_object in gvcfs) {
		File gvcf = gvcf_object.data
		File gvcf_index = gvcf_object.data_index
	}

	call PbsvCall.pbsv_call {
		input:
			sample_id = cohort_id,
			svsigs = svsigs,
			sample_count = sample_count,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			pbsv_call_mem_gb = pbsv_call_mem_gb,
			runtime_attributes = default_runtime_attributes
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = pbsv_call.pbsv_vcf,
			runtime_attributes = default_runtime_attributes
	}

	call Glnexus.glnexus {
		input:
			cohort_id = cohort_id,
			gvcfs = gvcf,
			gvcf_indices = gvcf_index,
			reference_name = reference.name,
			glnexus_mem_gb = glnexus_mem_gb,
			runtime_attributes = default_runtime_attributes
	}

	call PhaseVcf.phase_vcf {
		input:
			vcf = {"data": glnexus.vcf, "data_index": glnexus.vcf_index},
			aligned_bams = aligned_bams,
			reference_fasta = reference.fasta,
			reference_chromosome_lengths = reference.chromosome_lengths,
			regions = reference.chromosomes,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData sv_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		IndexData phased_joint_called_vcf = phase_vcf.phased_vcf
		File whatshap_stats_gtf = phase_vcf.whatshap_stats_gtf
		File whatshap_stats_tsv = phase_vcf.whatshap_stats_tsv
		File whatshap_stats_blocklist = phase_vcf.whatshap_stats_blocklist
	}

	parameter_meta {
		cohort_id: {help: "Sample ID for the cohort; used for naming files"}
		aligned_bams: {help: "Bam and index aligned to the reference genome for each movie associated with all samples in the cohort"}
		svsigs: {help: "pbsv svsig files for each sample and movie bam in the cohort"}
		gvcfs: {help: "gVCF for each sample in the cohort"}
		reference: {help: "Reference genome data"}
		pbsv_call_mem_gb: {help: "Optional amount of RAM in GB for pbsv_call"}
		glnexus_mem_gb: {help: "Optional amount of RAM in GB for glnexus"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}
