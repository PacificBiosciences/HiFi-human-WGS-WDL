version 1.0

# Run joint genotyping for a cohort. This workflow will be run if there is more than one sample in the cohort.

import "../humanwgs_structs.wdl"
import "../wdl-common/wdl/tasks/pbsv_call.wdl" as PbsvCall
import "../wdl-common/wdl/tasks/concat_vcf.wdl" as ConcatVcf
import "../wdl-common/wdl/tasks/glnexus.wdl" as Glnexus
import "../wdl-common/wdl/workflows/hiphase/hiphase.wdl" as HiPhase

workflow cohort_analysis {
	input {
		String cohort_id
		Array[String] sample_ids
		Array[IndexData] aligned_bams
		Array[File] svsigs
		Array[IndexData] gvcfs

		ReferenceData reference

		Int? pbsv_call_mem_gb
		Int? glnexus_mem_gb

		RuntimeAttributes default_runtime_attributes
	}

	Int sample_count = length(sample_ids)
	Array[Array[String]] pbsv_splits = read_json(reference.pbsv_splits)

	scatter (gvcf_object in gvcfs) {
		File gvcf = gvcf_object.data
		File gvcf_index = gvcf_object.data_index
	}

	scatter (shard_index in range(length(pbsv_splits))) {
        Array[String] region_set = pbsv_splits[shard_index]

		call PbsvCall.pbsv_call {
			input:
				sample_id = cohort_id + ".joint",
				svsigs = svsigs,
				sample_count = sample_count,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				reference_name = reference.name,
				shard_index = shard_index,
				regions = region_set,
				mem_gb = pbsv_call_mem_gb,
				runtime_attributes = default_runtime_attributes
		}
	}

	call ConcatVcf.concat_vcf {
		input:
			vcfs = pbsv_call.pbsv_vcf,
			vcf_indices = pbsv_call.pbsv_vcf_index,
			output_vcf_name = "~{cohort_id}.joint.~{reference.name}.pbsv.vcf.gz",
			runtime_attributes = default_runtime_attributes
	}

	IndexData zipped_pbsv_vcf = {
		"data": concat_vcf.concatenated_vcf,
		"data_index": concat_vcf.concatenated_vcf_index
	}

	call Glnexus.glnexus {
		input:
			cohort_id = cohort_id + ".joint",
			gvcfs = gvcf,
			gvcf_indices = gvcf_index,
			reference_name = reference.name,
			mem_gb = glnexus_mem_gb,
			runtime_attributes = default_runtime_attributes
	}

	IndexData glnexus_vcf = {
		"data": glnexus.vcf,
		"data_index": glnexus.vcf_index
	}

	call HiPhase.hiphase {
		# VCF order: small variants, SVs
		input:
			id = cohort_id + ".joint",
			refname = reference.name,
			sample_ids = sample_ids,
			vcfs = [glnexus_vcf, zipped_pbsv_vcf],
			bams = aligned_bams,
			haplotag = false,
			reference_fasta = reference.fasta,
			default_runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData phased_joint_small_variant_vcf = hiphase.phased_vcfs[0]
		IndexData phased_joint_sv_vcf = hiphase.phased_vcfs[1]
		File hiphase_stats = hiphase.hiphase_stats
		File hiphase_blocks = hiphase.hiphase_blocks
	}

	parameter_meta {
		cohort_id: {help: "Cohort ID; used for naming files"}
		sample_ids: {help: "Sample IDs for all samples in the cohort"}
		aligned_bams: {help: "BAM and index aligned to the reference genome for each movie associated with all samples in the cohort"}
		svsigs: {help: "pbsv svsig files for each sample and movie BAM in the cohort"}
		gvcfs: {help: "gVCF for each sample in the cohort"}
		reference: {help: "Reference genome data"}
		pbsv_call_mem_gb: {help: "Optional amount of RAM in GB for pbsv_call; default 64 for cohorts N<=3, 96 for cohorts N>3"}
		glnexus_mem_gb: {help: "Optional amount of RAM in GB for glnexus; default 30"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}
