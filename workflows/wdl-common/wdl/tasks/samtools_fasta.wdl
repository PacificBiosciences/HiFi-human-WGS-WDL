version 1.0

# Convert a BAM to a FASTA file using samtools

import "../structs.wdl"

task samtools_fasta {
	input {
		File bam

		RuntimeAttributes runtime_attributes
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 4
	Int disk_size = ceil(size(bam, "GB") * 3.5 + 20)

	command <<<
		set -euo pipefail

		samtools --version

		samtools fasta \
			-@ ~{threads - 1} \
			~{bam} \
		> ~{bam_basename}.fasta
	>>>

	output {
		File reads_fasta = "~{bam_basename}.fasta"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:83ca955c4a83f72f2cc229f41450eea00e73333686f3ed76f9f4984a985c85bb"
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
