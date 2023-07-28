version 1.0

# Run joint calling using GLnexus

import "../structs.wdl"

task glnexus {
	input {
		String cohort_id
		Array[File] gvcfs
		Array[File] gvcf_indices

		String reference_name

		File? regions_bed

		Int mem_gb = 30

		RuntimeAttributes runtime_attributes
	}

	Int threads = 24
	Int disk_size = ceil((size(gvcfs[0], "GB") * length(gvcfs)) * 2 + 100)

	command <<<
		set -euo pipefail

		# glneux_cli has no version option
		glnexus_cli --help 2>&1 | grep -Eo 'glnexus_cli release v[0-9a-f.-]+'

		glnexus_cli \
			--threads ~{threads} \
			--mem-gbytes ~{mem_gb} \
			--dir ~{cohort_id}.~{reference_name}.GLnexus.DB \
			--config DeepVariant_unfiltered \
			~{"--bed " + regions_bed} \
			~{sep=' ' gvcfs} \
		> ~{cohort_id}.~{reference_name}.deepvariant.glnexus.bcf

		bcftools --version

		bcftools view \
			--threads ~{threads} \
			--output-type z \
			--output-file ~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz \
			~{cohort_id}.~{reference_name}.deepvariant.glnexus.bcf

		tabix --version

		tabix ~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz
	>>>

	output {
		File vcf = "~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz"
		File vcf_index = "~{cohort_id}.~{reference_name}.deepvariant.glnexus.vcf.gz.tbi"
	}

	runtime {
		docker: "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
