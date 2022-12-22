version 1.0

import "../common/structs.wdl"

workflow phase_vcf {
	input {
		IndexData vcf
		Array[IndexData] aligned_bams

		ReferenceData reference

		String container_registry
		Boolean preemptible
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
				container_registry = container_registry,
				preemptible = preemptible
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
				container_registry = container_registry,
				preemptible = preemptible
		}
	}

	call bcftools_concat {
		input:
			vcfs = whatshap_phase.phased_vcf,
			vcf_indices = whatshap_phase.phased_vcf_index,
			output_vcf_name = "~{vcf_basename}.phased.vcf.gz",
			container_registry = container_registry,
			preemptible = preemptible
	}

	call whatshap_stats {
		input:
			phased_vcf = bcftools_concat.concatenated_vcf,
			phased_vcf_index = bcftools_concat.concatenated_vcf_index,
			reference_chromosome_lengths = reference.chromosome_lengths,
			container_registry = container_registry,
			preemptible = preemptible
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
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task split_vcf {
	input {
		File vcf
		File vcf_index
		String region

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/htslib:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
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

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/whatshap:b1a46c6"
		cpu: 1
		memory: "8 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task bcftools_concat {
	input {
		Array[File] vcfs
		Array[File] vcf_indices
		String output_vcf_name

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task whatshap_stats {
	input {
		File phased_vcf
		File phased_vcf_index

		File reference_chromosome_lengths

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/whatshap:b1a46c6"
		cpu: 1
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
