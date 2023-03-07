version 1.0

import "../common/structs.wdl"

workflow phase_vcf {
	input {
		IndexData vcf
		Array[IndexData] aligned_bams

		ReferenceData reference

		RuntimeAttributes default_runtime_attributes
	}

	String vcf_basename = basename(vcf.data, ".vcf.gz")

	scatter (bam_object in aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.data_index
	}

	scatter (chromosome in reference.chromosomes) {
		call split_vcf {
			input:
				vcf = vcf.data,
				vcf_index = vcf.data_index,
				region = chromosome,
				runtime_attributes = default_runtime_attributes
		}

		call whatshap_phase {
			input:
				vcf = split_vcf.region_vcf,
				vcf_index = split_vcf.region_vcf_index,
				chromosome = chromosome,
				aligned_bams = aligned_bam,
				aligned_bam_indices = aligned_bam_index,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				runtime_attributes = default_runtime_attributes
		}
	}

	call bcftools_concat {
		input:
			vcfs = whatshap_phase.phased_vcf,
			vcf_indices = whatshap_phase.phased_vcf_index,
			output_vcf_name = "~{vcf_basename}.phased.vcf.gz",
			runtime_attributes = default_runtime_attributes
	}

	call whatshap_stats {
		input:
			phased_vcf = bcftools_concat.concatenated_vcf,
			phased_vcf_index = bcftools_concat.concatenated_vcf_index,
			reference_chromosome_lengths = reference.chromosome_lengths,
			runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData phased_vcf = {"data": bcftools_concat.concatenated_vcf, "data_index": bcftools_concat.concatenated_vcf_index}
		File whatshap_stats_gtf = whatshap_stats.gtf
		File whatshap_stats_tsv = whatshap_stats.tsv
		File whatshap_stats_blocklist = whatshap_stats.blocklist
	}

	parameter_meta {
		vcf: {help: "VCF to phase"}
		aligned_bams: {help: "Bam and index aligned to the reference genome for each movie associated with the sample"}
		reference: {help: "Reference genome data"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task split_vcf {
	input {
		File vcf
		File vcf_index
		String region

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		tabix \
			-h \
			~{vcf} \
			~{region} \
		> ~{vcf_basename}.~{region}.vcf

		bgzip ~{vcf_basename}.~{region}.vcf
		tabix ~{vcf_basename}.~{region}.vcf.gz
	>>>

	output {
		File region_vcf = "~{vcf_basename}.~{region}.vcf.gz"
		File region_vcf_index = "~{vcf_basename}.~{region}.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/htslib:1.14"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task whatshap_phase {
	input {
		File vcf
		File vcf_index
		String chromosome

		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference
		File reference_index

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(aligned_bams[0], "GB") * length(aligned_bams)) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap phase \
			--indels \
			--reference ~{reference} \
			--chromosome ~{chromosome} \
			--output ~{vcf_basename}.phased.vcf.gz \
			~{vcf} \
			~{sep=' ' aligned_bams}

		tabix ~{vcf_basename}.phased.vcf.gz
	>>>

	output {
		File phased_vcf = "~{vcf_basename}.phased.vcf.gz"
		File phased_vcf_index = "~{vcf_basename}.phased.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/whatshap:1.4"
		cpu: 1
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

task bcftools_concat {
	input {
		Array[File] vcfs
		Array[File] vcf_indices
		String output_vcf_name

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil(size(vcfs[0], "GB") * length(vcfs) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools concat \
			--allow-overlaps \
			--output ~{output_vcf_name} \
			--output-type z \
			~{sep=' ' vcfs}

		tabix "~{output_vcf_name}"
	>>>

	output {
		File concatenated_vcf = "~{output_vcf_name}"
		File concatenated_vcf_index = "~{output_vcf_name}.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools:1.14"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task whatshap_stats {
	input {
		File phased_vcf
		File phased_vcf_index

		File reference_chromosome_lengths

		RuntimeAttributes runtime_attributes
	}

	String output_basename = basename(phased_vcf, ".vcf.gz")
	Int disk_size = ceil(size(phased_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap stats \
			--gtf ~{output_basename}.gtf \
			--tsv ~{output_basename}.tsv \
			--block-list ~{output_basename}.blocklist \
			--chr-lengths ~{reference_chromosome_lengths} \
			~{phased_vcf}
	>>>

	output {
		File gtf = "~{output_basename}.gtf"
		File tsv = "~{output_basename}.tsv"
		File blocklist = "~{output_basename}.blocklist"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/whatshap:1.4"
		cpu: 1
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
