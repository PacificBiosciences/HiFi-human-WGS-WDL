version 1.0

import "../common/structs.wdl"
import "../common/tasks/pbsv_call.wdl" as PbsvCall
import "../common/tasks/zip_index_vcf.wdl" as ZipIndexVcf
import "../phase_vcf/phase_vcf.wdl" as PhaseVcf

workflow cohort_analysis {
	input {
		String cohort_id
		Array[IndexData] aligned_bams
		Array[File] svsigs
		Array[IndexData] gvcfs

		ReferenceData reference

		RuntimeAttributes spot_runtime_attributes
	}

	scatter (gvcf_object in gvcfs) {
		File gvcf = gvcf_object.data
		File gvcf_index = gvcf_object.data_index
	}

	call PbsvCall.pbsv_call {
		input:
			sample_id = cohort_id,
			svsigs = svsigs,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			runtime_attributes = spot_runtime_attributes
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = pbsv_call.pbsv_vcf,
			runtime_attributes = spot_runtime_attributes
	}

	call glnexus {
		input:
			cohort_id = cohort_id,
			gvcfs = gvcf,
			gvcf_indices = gvcf_index,
			reference_name = reference.name,
			runtime_attributes = spot_runtime_attributes
	}

	call PhaseVcf.phase_vcf {
		input:
			vcf = {"data": glnexus.vcf, "data_index": glnexus.vcf_index},
			aligned_bams = aligned_bams,
			reference = reference,
			spot_runtime_attributes = spot_runtime_attributes
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
		spot_runtime_attributes: {help: "RuntimeAttributes for spot (preemptible) tasks"}
	}
}

task glnexus {
	input {
		String cohort_id
		Array[File] gvcfs
		Array[File] gvcf_indices

		String reference_name

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int mem_gbytes = 30
	Int disk_size = ceil((size(gvcfs[0], "GB") * length(gvcfs)) * 2 + 100)

	command <<<
		set -euo pipefail

		glnexus_cli \
			--threads ~{threads} \
			--mem-gbytes ~{mem_gbytes} \
			--dir ~{cohort_id}.~{reference_name}.GLnexus.DB \
			--config DeepVariant_unfiltered \
			~{sep=' ' gvcfs} \
		> ~{cohort_id}.~{reference_name}.deepvariant.glnexus.bcf

		bcftools view \
			--threads ~{threads} \
			--output-type z \
			--output-file ~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz \
			~{cohort_id}.~{reference_name}.deepvariant.glnexus.bcf

		tabix ~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz
	>>>

	output {
		File vcf = "~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz"
		File vcf_index = "~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz.tbi"
	}

	runtime {
		docker: "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
		cpu: threads
		memory: mem_gbytes + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
