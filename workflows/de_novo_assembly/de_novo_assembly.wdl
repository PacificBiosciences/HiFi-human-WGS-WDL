version 1.0

import "../common/common.wdl" as common
import "../common/structs.wdl"

workflow de_novo_assembly {
	input {
		Sample sample

		ReferenceData reference

		String container_registry
	}

	scatter (movie_bam in sample.movie_bams) {
		call samtools_fasta {
			input:
				bam = movie_bam.data,
				container_registry = container_registry
		}
	}

	call hifiasm_assemble {
		input:
			sample_id = sample.sample_id,
			reads_fastas = samtools_fasta.reads_fasta,
			container_registry = container_registry
	}

	scatter (gfa in hifiasm_assemble.assembly_hap_gfas) {
		call gfa2fa {
			input:
				gfa = gfa,
				container_registry = container_registry
		}

		call bgzip_fasta {
			input:
				fasta = gfa2fa.fasta,
				container_registry = container_registry
		}

		call asm_stats {
			input:
				zipped_fasta = bgzip_fasta.zipped_fasta,
				reference_index = reference.fasta.data_index,
				container_registry = container_registry
		}
	}

	call align_hifiasm {
		input:
			sample_id = sample.sample_id,
			query_sequences = bgzip_fasta.zipped_fasta,
			reference = reference.fasta.data,
			reference_name = reference.name,
			container_registry = container_registry
	}

	call htsbox {
		input:
			bam = align_hifiasm.asm_bam,
			bam_index = align_hifiasm.asm_bam_index,
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	call common.zip_index_vcf {
		input:
			vcf = htsbox.htsbox_vcf,
			container_registry = container_registry
	}

	call bcftools_stats {
		input:
			vcf = zip_index_vcf.zipped_vcf,
			bam = align_hifiasm.asm_bam,
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	output {
		Array[File] assembly_noseq_gfas = hifiasm_assemble.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = hifiasm_assemble.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = bgzip_fasta.zipped_fasta
		Array[File] assembly_stats = asm_stats.assembly_stats
		IndexData asm_bam = {"data": align_hifiasm.asm_bam, "data_index": align_hifiasm.asm_bam_index}
		IndexData htsbox_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		File htsbox_vcf_stats = bcftools_stats.vcf_stats
	}

	parameter_meta {
		sample: {help: "Sample ID and unaligned movie bams and indices associated with the sample"}
		reference: {help: "ReferenceData"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

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

task hifiasm_assemble {
	input {
		String sample_id
		Array[File] reads_fastas

		String container_registry
	}

	String prefix = "~{sample_id}.asm"
	Int threads = 48
	Int disk_size = ceil((size(reads_fastas[0], "GB") * length(reads_fastas)) * 2 + 20)

	command <<<
		set -euo pipefail

		hifiasm \
			-o ~{prefix} \
			-t ~{threads} \
			~{sep=' ' reads_fastas}
	>>>

	output {
		Array[File] assembly_hap_gfas = glob("~{prefix}.bp.hap[12].p_ctg.gfa")
		Array[File] assembly_noseq_gfas = glob("~{prefix}.bp.hap[12].p_ctg.noseq.gfa")
		Array[File] assembly_lowQ_beds = glob("~{prefix}.bp.hap[12].p_ctg.lowQ.bed")
	}

	runtime {
		docker: "~{container_registry}/hifiasm:b1a46c6"
		cpu: threads
		memory: "192 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task gfa2fa {
	input {
		File gfa

		String container_registry
	}

	String gfa_basename = basename(gfa, ".gfa")
	Int disk_size = ceil(size(gfa, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		gfatools gfa2fa \
			~{gfa} \
		> ~{gfa_basename}.fasta
	>>>

	output {
		File fasta = "~{gfa_basename}.fasta"
	}

	runtime {
		docker: "~{container_registry}/gfatools:b1a46c6"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bgzip_fasta {
	input {
		File fasta

		String container_registry
	}

	String fasta_basename = basename(fasta)
	Int threads = 4
	Int disk_size = ceil(size(fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		bgzip \
			--threads ~{threads} \
			--stdout \
			~{fasta} \
		> ~{fasta_basename}.gz
	>>>

	output {
		File zipped_fasta = "~{fasta_basename}.gz"
	}

	runtime {
		docker: "~{container_registry}/htslib:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task asm_stats {
	input {
		File zipped_fasta

		File reference_index

		String container_registry
	}

	String zipped_fasta_basename = basename(zipped_fasta, ".gz")
	Int disk_size = ceil(size(zipped_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		k8 \
			/opt/scripts/calN50/calN50.js \
			-f ~{reference_index} \
			~{zipped_fasta} \
		> ~{zipped_fasta_basename}.stats.txt
	>>>

	output {
		File assembly_stats = "~{zipped_fasta_basename}.stats.txt"
	}

	runtime {
		docker: "~{container_registry}/k8:b1a46c6"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task align_hifiasm {
	input {
		String sample_id
		Array[File] query_sequences

		File reference
		String reference_name

		String container_registry
	}

	Int threads = 16
	Int disk_size = ceil((size(query_sequences[0], "GB") * length(query_sequences) + size(reference, "GB")) * 2 + 20)

	command <<<
		minimap2 \
			-t ~{threads - 4} \
			-L \
			--secondary=no \
			--eqx \
			-a \
			-x asm5 \
			-R "@RG\\tID:~{sample_id}_hifiasm\\tSM:~{sample_id}" \
			~{sep=' ' query_sequences} \
		| samtools sort \
			-@ 4 \
			-T ./TMP \
			-m 8G \
			-O BAM \
			-o ~{sample_id}.asm.~{reference_name}.bam

		samtools index ~{sample_id}.asm.~{reference_name}.bam
	>>>

	output {
		File asm_bam = "~{sample_id}.asm.~{reference_name}.bam"
		File asm_bam_index = "~{sample_id}.asm.~{reference_name}.bam.bai"
	}

	runtime {
		docker: "~{container_registry}/align_hifiasm:b1a46c6"
		cpu: threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
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
