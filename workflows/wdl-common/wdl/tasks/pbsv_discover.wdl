version 1.0

# Generate svsigs to be used in SV calling using pbsv

import "../structs.wdl"

task pbsv_discover {
	input {
		File aligned_bam
		File aligned_bam_index

		File reference_tandem_repeat_bed

		RuntimeAttributes runtime_attributes
	}

	String prefix = basename(aligned_bam, ".bam")
	Int disk_size = ceil((size(aligned_bam, "GB") + size(reference_tandem_repeat_bed, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv --version

		pbsv discover \
			--log-level INFO \
			--hifi \
			--tandem-repeats ~{reference_tandem_repeat_bed} \
			~{aligned_bam} \
			~{prefix}.svsig.gz
	>>>

	output {
		File svsig = "~{prefix}.svsig.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pbsv@sha256:d78ee6deb92949bdfde98d3e48dab1d871c177d48d8c87c73d12c45bdda43446"
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
