version 1.0

task pbsv_call {
	input {
		String sample_id
		Array[File] svsigs

		File reference
		File reference_index
		String reference_name

		String container_registry
		Boolean preemptible
	}

	# TODO does this slow it down? 2:15 at 8 cores
	Int threads = 4
	Int disk_size = ceil((size(svsigs[0], "GB") * length(svsigs) + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv call \
			--hifi \
			--min-sv-length 20 \
			--num-threads ~{threads} \
			~{reference} \
			~{sep=' ' svsigs} \
			~{sample_id}.~{reference_name}.pbsv.vcf
	>>>

	output {
		File pbsv_vcf = "~{sample_id}.~{reference_name}.pbsv.vcf"
	}

	runtime {
		docker: "~{container_registry}/pbsv:b1a46c6"
		cpu: threads
		memory: "56 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
