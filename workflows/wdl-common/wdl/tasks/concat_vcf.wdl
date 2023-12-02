version 1.0

# Concatenate and sort VCFs

import "../structs.wdl"

task concat_vcf {
	input {
		Array[File] vcfs
		Array[File] vcf_indices

		String output_vcf_name

		RuntimeAttributes runtime_attributes
	}

	Int threads = 4
	Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)
	
	command <<<
		set -euo pipefail

		mkdir vcfs
		while read -r input || [[ -n "${input}" ]]; do
			ln -s "${input}" vcfs
		done < ~{write_lines(flatten([vcfs,vcf_indices]))}

		find vcfs -name "*.vcf.gz" > vcf.list

		bcftools --version

		bcftools concat \
			--allow-overlaps \
			--threads ~{threads - 1} \
			--output-type z \
			--output ~{output_vcf_name} \
			--file-list vcf.list

		bcftools index --tbi ~{output_vcf_name}
	>>>

	output {
		File concatenated_vcf = "~{output_vcf_name}"
		File concatenated_vcf_index = "~{output_vcf_name}.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
		cpu: threads
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
