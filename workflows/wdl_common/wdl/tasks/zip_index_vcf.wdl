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
		docker: "~{runtime_attributes.container_registry}/htslib@sha256:112ec32ffd4fc0891e88270f5a8fa2b2820809dee5b3ef7f27fc12934f43c0c1"
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
