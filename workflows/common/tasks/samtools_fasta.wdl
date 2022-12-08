version 1.0

task samtools_fasta {
	input {
		File bam

		String container_registry
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 4
	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		samtools fasta \
			-@ ~{threads - 1} \
			~{bam} \
		> ~{bam_basename}.fasta
	>>>

	output {
		File reads_fasta = "~{bam_basename}.fasta"
	}

	runtime {
		docker: "~{container_registry}/samtools:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}
