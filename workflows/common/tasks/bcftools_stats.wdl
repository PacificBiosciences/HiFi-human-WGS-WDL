version 1.0

task bcftools_stats {
	input {
		File vcf
		String? params

		File? reference
		File? bam

		String container_registry
		Boolean preemptible
	}

	String vcf_basename = basename(vcf, ".gz")

	Int threads = 2
	Int reference_size = if (defined(reference)) then ceil(size(reference, "GB")) else 0
	Int bam_size = if (defined(bam)) then ceil(size(bam, "GB")) else 0
	Int disk_size = ceil((size(vcf, "GB") + reference_size + bam_size) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools stats \
			--threads ~{threads} \
			~{params} \
			~{"--fasta-ref " + reference} \
			~{"--samples " + bam}
			~{vcf} \
		> ~{vcf_basename}.stats.txt
	>>>

	output {
		File stats = "~{vcf_basename}.stats.txt"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: threads
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
