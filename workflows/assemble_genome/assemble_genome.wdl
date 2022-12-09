version 1.0

import "../common/structs.wdl"

workflow assemble_genome {
	input {
		String sample_id
		Array[File] reads_fastas

		ReferenceData reference

		String? hifiasm_extra_params
		File? father_yak
		File? mother_yak

		String container_registry
	}

	call hifiasm_assemble {
		input:
			sample_id = sample_id,
			reads_fastas = reads_fastas,
			extra_params = hifiasm_extra_params,
			father_yak = father_yak,
			mother_yak = mother_yak,
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
			sample_id = sample_id,
			query_sequences = bgzip_fasta.zipped_fasta,
			reference = reference.fasta.data,
			reference_name = reference.name,
			container_registry = container_registry
	}

	output {
		Array[File] assembly_noseq_gfas = hifiasm_assemble.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = hifiasm_assemble.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = bgzip_fasta.zipped_fasta
		Array[File] assembly_stats = asm_stats.assembly_stats
		IndexData asm_bam = {"data": align_hifiasm.asm_bam, "data_index": align_hifiasm.asm_bam_index}
	}

	parameter_meta {
		sample_id: {help: "Sample ID; used for naming files"}
		reads_fastas: {help: "Reads in fasta format to be used for assembly; one for each movie bam to be used in assembly. Reads fastas from one or more sample may be combined to use in the assembly"}
		reference: {help: "Reference genome data"}
		hiiasm_extra_params: {help: "[OPTIONAL] Additional parameters to pass to hifiasm assembly"}
		father_yak: {help: "[OPTIONAL] kmer counts for the father; required if running trio-based assembly"}
		mother_yak: {help: "[OPTIONAL] kmer counts for the mother; required if running trio-based assembly"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task hifiasm_assemble {
	input {
		String sample_id
		Array[File] reads_fastas

		String? extra_params
		File? father_yak
		File? mother_yak

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
			~{extra_params} \
			~{"-1 " + father_yak} \
			~{"-2 " + mother_yak} \
			~{sep=' ' reads_fastas}
	>>>

	output {
		Array[File] assembly_hap_gfas = glob("~{prefix}.*.hap[12].p_ctg.gfa")
		Array[File] assembly_noseq_gfas = flatten([
			glob("~{prefix}.*.hap[12].p_ctg.noseq.gfa"),
			glob("~{prefix}.dip.[pr]_utg.noseq.gfa")
		])
		Array[File] assembly_lowQ_beds = flatten([
			glob("~{prefix}.*.hap[12].p_ctg.lowQ.bed"),
			glob("~{prefix}.dip.[pr]_utg.lowQ.bed")
		])
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
		set -euo pipefail

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
