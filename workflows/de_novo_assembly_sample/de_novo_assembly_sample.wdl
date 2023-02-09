version 1.0

import "../common/structs.wdl"
import "../common/tasks/samtools_fasta.wdl" as SamtoolsFasta
import "../assemble_genome/assemble_genome.wdl" as AssembleGenome
import "../common/tasks/zip_index_vcf.wdl" as ZipIndexVcf
import "../common/tasks/bcftools_stats.wdl" as BcftoolsStats

workflow de_novo_assembly_sample {
	input {
		Sample sample

		ReferenceData reference

		String container_registry
		Boolean preemptible
	}

	scatter (movie_bam in sample.movie_bams) {
		call SamtoolsFasta.samtools_fasta {
			input:
				bam = movie_bam,
				container_registry = container_registry,
				preemptible = preemptible
		}
	}

	call AssembleGenome.assemble_genome {
		input:
			sample_id = sample.sample_id,
			reads_fastas = samtools_fasta.reads_fasta,
			reference = reference,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call htsbox {
		input:
			bam = assemble_genome.asm_bam.data,
			bam_index = assemble_genome.asm_bam.data_index,
			reference = reference.fasta.data,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = htsbox.htsbox_vcf,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call BcftoolsStats.bcftools_stats {
		input:
			vcf = zip_index_vcf.zipped_vcf,
			params = "--samples ~{basename(assemble_genome.asm_bam.data)}",
			reference = reference.fasta.data,
			container_registry = container_registry,
			preemptible = preemptible
	}

	output {
		Array[File] assembly_noseq_gfas = assemble_genome.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = assemble_genome.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = assemble_genome.zipped_assembly_fastas
		Array[File] assembly_stats = assemble_genome.assembly_stats
		IndexData asm_bam = assemble_genome.asm_bam
		IndexData htsbox_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		File htsbox_vcf_stats = bcftools_stats.stats
	}

	parameter_meta {
		sample: {help: "Sample information and associated data files"}
		reference: {help: "Reference genome data"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task htsbox {
	input {
		File bam
		File bam_index

		File reference

		String container_registry
		Boolean preemptible
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 2
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 3 + 200)

	command <<<
		set -euo pipefail

		# Ensure the sample is named based on the bam basename (not the full path)
		cp ~{bam} .

		htsbox pileup \
			-q20 \
			-c \
			-f ~{reference} \
			~{basename(bam)} \
		> ~{bam_basename}.htsbox.vcf
	>>>

	output {
		File htsbox_vcf = "~{bam_basename}.htsbox.vcf"
	}

	runtime {
		docker: "~{container_registry}/htsbox:r346"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
