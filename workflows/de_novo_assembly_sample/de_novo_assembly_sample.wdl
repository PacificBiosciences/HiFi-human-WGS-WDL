version 1.0

import "../common/structs.wdl"
import "../common/tasks/samtools_fasta.wdl" as SamtoolsFasta
import "../assemble_genome/assemble_genome.wdl" as AssembleGenome
import "../common/tasks/zip_index_vcf.wdl" as ZipIndexVcf

workflow de_novo_assembly_sample {
	input {
		Sample sample

		ReferenceData reference

		String container_registry
	}

	scatter (movie_bam in sample.movie_bams) {
		call SamtoolsFasta.samtools_fasta {
			input:
				bam = movie_bam,
				container_registry = container_registry
		}
	}

	call AssembleGenome.assemble_genome {
		input:
			sample_id = sample.sample_id,
			reads_fastas = samtools_fasta.reads_fasta,
			reference = reference,
			container_registry = container_registry
	}

	call htsbox {
		input:
			bam = assemble_genome.asm_bam.data,
			bam_index = assemble_genome.asm_bam.data_index,
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = htsbox.htsbox_vcf,
			container_registry = container_registry
	}

	call bcftools_stats {
		input:
			vcf = zip_index_vcf.zipped_vcf,
			bam = assemble_genome.asm_bam.data,
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	output {
		Array[File] assembly_noseq_gfas = assemble_genome.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = assemble_genome.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = assemble_genome.zipped_assembly_fastas
		Array[File] assembly_stats = assemble_genome.assembly_stats
		IndexData asm_bam = assemble_genome.asm_bam
		IndexData htsbox_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		File htsbox_vcf_stats = bcftools_stats.vcf_stats
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
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 4
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		htsbox pileup \
			-q20 \
			-c \
			-f \
			~{reference} \
			~{bam} \
		> ~{bam_basename}.htsbox.vcf
	>>>

	output {
		File htsbox_vcf = "~{bam_basename}.htsbox.vcf"
	}

	runtime {
		docker: "~{container_registry}/htsbox:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcftools_stats {
	input {
		File vcf
		File bam

		File reference

		String container_registry
	}

	String vcf_basename = basename(vcf, ".gz")
	Int threads = 4
	Int disk_size = ceil((size(vcf, "GB") + size(bam, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools stats \
			--threads ~{threads - 1} \
			--fasta-ref ~{reference} \
			--samples ~{bam} \
			~{vcf} \
		> ~{vcf_basename}.stats.txt
	>>>

	output {
		File vcf_stats = "~{vcf_basename}.stats.txt"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}