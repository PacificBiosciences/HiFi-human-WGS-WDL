version 1.0

task mosdepth {
	input {
		File aligned_bam
		File aligned_bam_index

		String container_registry
		Boolean preemptible
	}

	String prefix = basename(aligned_bam, ".bam")
	# TODO does this slow down at 2 cores instead of 4? 6 mins at 4
	Int threads = 2
	Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		mosdepth \
			--threads ~{threads} \
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
		docker: "~{container_registry}/mosdepth:b1a46c6"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
