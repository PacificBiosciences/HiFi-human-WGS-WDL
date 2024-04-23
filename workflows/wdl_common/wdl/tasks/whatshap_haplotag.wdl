version 1.0

# Haplotag an aligned BAM file using a phased VCF with WhatsHap

import "../structs.wdl"

task whatshap_haplotag {
	input {
		File phased_vcf
		File phased_vcf_index

		File aligned_bam
		File aligned_bam_index

		File reference
		File reference_index

		String? params
		String? output_bam_name

		RuntimeAttributes runtime_attributes
	}

	String output_bam = select_first([output_bam_name, "~{basename(aligned_bam, '.bam')}.haplotagged.bam"])
	Int threads = 4
	Int disk_size = ceil((size(phased_vcf, "GB") + size(aligned_bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap --version

		whatshap haplotag \
			~{params} \
			--tag-supplementary \
			--output-threads ~{threads} \
			--reference ~{reference} \
			--output ~{output_bam} \
			~{phased_vcf} \
			~{aligned_bam}

		samtools --version

		samtools index \
			-@ ~{threads - 1} \
			~{output_bam}
	>>>

	output {
		File haplotagged_bam = "~{output_bam}"
		File haplotagged_bam_index = "~{output_bam}.bai"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/whatshap@sha256:cf889352d7b5b977e70bd6386f36017006c093d5f6fffca4c8ee8eeda4d10b95"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
