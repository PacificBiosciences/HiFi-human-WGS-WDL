version 1.0

task zip_index_vcf {
	input {
		File vcf

		String container_registry
		Boolean preemptible
	}

	String vcf_basename = basename(vcf)
	Int threads = 4
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		bgzip \
			--threads ~{threads} \
			--stdout \
			~{vcf} \
		> ~{vcf_basename}.gz

		tabix \
			--preset vcf \
			~{vcf_basename}.gz
	>>>

	output {
		File zipped_vcf = "~{vcf_basename}.gz"
		File zipped_vcf_index = "~{vcf_basename}.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/htslib:1.14"
		cpu: threads
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
