version 1.0

# Calculate VCF stats

import "../structs.wdl"

task bcftools_stats {
	input {
		File vcf
		String? params

		File? reference

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf, ".gz")

	Int threads = 2
	Int reference_size = if (defined(reference)) then ceil(size(reference, "GB")) else 0
	Int disk_size = ceil((size(vcf, "GB") + reference_size) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools --help

		bcftools stats \
			--threads ~{threads - 1} \
			~{params} \
			~{"--fasta-ref " + reference} \
			~{vcf} \
		> ~{vcf_basename}.stats.txt
	>>>

	output {
		File stats = "~{vcf_basename}.stats.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
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
