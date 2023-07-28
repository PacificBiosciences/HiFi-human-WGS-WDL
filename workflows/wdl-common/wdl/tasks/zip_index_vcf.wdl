version 1.0

# Zip and index a VCF file

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

		bgzip --version

		bgzip \
			--threads ~{threads} \
			--stdout \
			~{vcf} \
		> ~{vcf_basename}.gz

		tabix --version

		tabix \
			--preset vcf \
			~{vcf_basename}.gz
	>>>

	output {
		File zipped_vcf = "~{vcf_basename}.gz"
		File zipped_vcf_index = "~{vcf_basename}.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/htslib@sha256:24ae834b9d4ba3ea3c23d77b2ce49b3a56a6e32d1367470e8e1160eb645019a9"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " LOCAL"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
