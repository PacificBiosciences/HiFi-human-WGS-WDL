version 1.0

# Phase and calculate stats for a VCF using WhatsHap.

import "../../structs.wdl"
import "../../tasks/whatshap_phase.wdl" as WhatshapPhase
import "../../tasks/whatshap_stats.wdl" as WhatshapStats

workflow phase_vcf {
	input {
		IndexData vcf
		Array[IndexData] aligned_bams

		IndexData reference_fasta
		File reference_chromosome_lengths
		Array[String] regions

		RuntimeAttributes default_runtime_attributes
	}

	String vcf_basename = basename(vcf.data, ".vcf.gz")

	scatter (bam_object in aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.data_index
	}

	scatter (region in regions) {
		call split_vcf {
			input:
				vcf = vcf.data,
				vcf_index = vcf.data_index,
				region = region,
				runtime_attributes = default_runtime_attributes
		}

		String chromosome = sub(region, ":.*", "")

		call WhatshapPhase.whatshap_phase {
			input:
				vcf = split_vcf.region_vcf,
				vcf_index = split_vcf.region_vcf_index,
				chromosome = chromosome,
				aligned_bams = aligned_bam,
				aligned_bam_indices = aligned_bam_index,
				reference = reference_fasta.data,
				reference_index = reference_fasta.data_index,
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

	call WhatshapStats.whatshap_stats {
		input:
			phased_vcf = bcftools_concat.concatenated_vcf,
			phased_vcf_index = bcftools_concat.concatenated_vcf_index,
			reference_chromosome_lengths = reference_chromosome_lengths,
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
		reference: {help: "Reference genome fasta and index"}
		reference_chromosome_lengths: {help: "File specifying the lengths of each of the reference chromosomes"}
		regions: {help: "Array of regions to run phasing on; can be in the format chr or chr:start-stop"}
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
	String region_substituted = sub(region, ":", "_")
	Int threads = 2
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		tabix --version

		tabix \
			-h \
			~{vcf} \
			~{region} \
		> ~{vcf_basename}.~{region_substituted}.vcf

		bgzip --version

		bgzip -@{threads} ~{vcf_basename}.~{region_substituted}.vcf
		tabix ~{vcf_basename}.~{region_substituted}.vcf.gz
	>>>

	output {
		File region_vcf = "~{vcf_basename}.~{region_substituted}.vcf.gz"
		File region_vcf_index = "~{vcf_basename}.~{region_substituted}.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/htslib@sha256:112ec32ffd4fc0891e88270f5a8fa2b2820809dee5b3ef7f27fc12934f43c0c1"
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

task bcftools_concat {
	input {
		Array[File] vcfs
		Array[File] vcf_indices
		String output_vcf_name

		RuntimeAttributes runtime_attributes
	}

	Int threads = 2
	Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools --version
		
		tabix --version

		bcftools concat \
			--threads ~{threads - 1} \
			--no-version \
			--file-list ~{write_lines(vcfs)} \
			--allow-overlaps \
			--output ~{output_vcf_name} \
			--output-type z

		tabix "~{output_vcf_name}"
	>>>

	output {
		File concatenated_vcf = "~{output_vcf_name}"
		File concatenated_vcf_index = "~{output_vcf_name}.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools@sha256:46720a7ab5feba5be06d5269454a6282deec13060e296f0bc441749f6f26fdec"
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
