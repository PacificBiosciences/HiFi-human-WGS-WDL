version 1.0

# Calculate summary stats using mosdepth

import "../structs.wdl"

task mosdepth {
	input {
		File aligned_bam
		File aligned_bam_index

		RuntimeAttributes runtime_attributes
	}

	String prefix = basename(aligned_bam, ".bam")
	Int threads = 4
	Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		mosdepth --version

		mosdepth \
			--threads ~{threads - 1} \
			--by 500 \
			--no-per-base \
			--use-median \
			~{prefix} \
			~{aligned_bam}
	>>>

	output {
		File summary = "~{prefix}.mosdepth.summary.txt"
		File region_bed = "~{prefix}.regions.bed.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:35d5e02facf4f38742e5cae9e5fdd3807c2b431dd8d881fd246b55e6d5f7f600"
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
