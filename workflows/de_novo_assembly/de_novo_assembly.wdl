version 1.0

import "../common/structs.wdl"

workflow de_novo_assembly {
	input {
		Sample sample

		IndexData reference_genome

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

	scatter (gfa in [hifiasm_assemble.hap1_gfa, hifiasm_assemble.hap2_gfa]) {
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
				reference_index = reference_genome.data_index,
				container_registry = container_registry
		}
	}

	output {
		Array[File] zipped_assembly_fastas = bgzip_fasta.zipped_fasta
		Array[File] assembly_stats = asm_stats.assembly_stats
	}

	parameter_meta {
		sample: {help: "Sample ID and unaligned movie bams and indices associated with the sample"}
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
		File hap1_gfa = "~{prefix}.bp.hap1.p_ctg.gfa"
		File hap2_gfa = "~{prefix}.bp.hap2.p_ctg.gfa"
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
