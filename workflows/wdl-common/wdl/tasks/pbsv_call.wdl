version 1.0

# Call SVs using pbsv

import "../structs.wdl"

task pbsv_call {
	input {
		String sample_id
		Array[File] svsigs
		Int? sample_count

		File reference
		File reference_index
		String reference_name

		Int mem_gb = if select_first([sample_count, 1]) > 3 then 96 else 64

		RuntimeAttributes runtime_attributes
	}

	Int threads = 8
	Int disk_size = ceil((size(svsigs[0], "GB") * length(svsigs) + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv --version

		pbsv call \
			--hifi \
			--min-sv-length 20 \
			--log-level INFO \
			--num-threads ~{threads} \
			~{reference} \
			~{sep=' ' svsigs} \
			~{sample_id}.~{reference_name}.pbsv.vcf
	>>>

	output {
		File pbsv_vcf = "~{sample_id}.~{reference_name}.pbsv.vcf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pbsv@sha256:798ca327f653c4e666b9f7c6a09260a762eea4e7e6864f490a87ed4106a53b98"
		cpu: threads
		memory: "~{mem_gb} GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
