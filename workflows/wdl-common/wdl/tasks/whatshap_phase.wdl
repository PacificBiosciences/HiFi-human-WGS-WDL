version 1.0

# Phase a VCF using WhatsHap

import "../structs.wdl"

task whatshap_phase {
	input {
		File vcf
		File vcf_index
		String? chromosome

		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference
		File reference_index

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(aligned_bams[0], "GB") * length(aligned_bams)) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap --version

		whatshap phase \
			--indels \
			--reference ~{reference} \
			~{"--chromosome " + chromosome} \
			--output ~{vcf_basename}.phased.vcf.gz \
			~{vcf} \
			~{sep=' ' aligned_bams}

		tabix ~{vcf_basename}.phased.vcf.gz
	>>>

	output {
		File phased_vcf = "~{vcf_basename}.phased.vcf.gz"
		File phased_vcf_index = "~{vcf_basename}.phased.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/whatshap@sha256:34957019d127e9c9c888a38061b28af8c1a42ec9e131bf1b806f70c6e96a1fca"
		cpu: 2
		memory: "8 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
