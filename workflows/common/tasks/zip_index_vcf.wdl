version 1.0

import "../structs.wdl"

task zip_index_vcf {
	input {
		File vcf

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf)
	Int threads = 4
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		bgzip \
			--threads ~{threads} \
			--stdout \
			~{vcf} \
		> ~{vcf_basename}.gz

		tabix \
			--preset vcf \
			~{vcf_basename}.gz
	>>>

	output {
		File zipped_vcf = "~{vcf_basename}.gz"
		File zipped_vcf_index = "~{vcf_basename}.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/htslib:b1a46c6"
		cpu: threads
		memory: "1 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
