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
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:cbe496e16773d4ad6f2eec4bd1b76ff142795d160f9dd418318f7162dcdaa685"
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
