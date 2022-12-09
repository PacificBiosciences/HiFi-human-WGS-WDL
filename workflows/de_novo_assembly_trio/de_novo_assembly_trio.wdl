version 1.0

import "../common/structs.wdl"
import "../common/tasks/samtools_fasta.wdl" as SamtoolsFasta
import "../assemble_genome/assemble_genome.wdl" as AssembleGenome

workflow de_novo_assembly_trio {
	input {
		Cohort cohort

		ReferenceData reference

		String container_registry
	}

	call parse_trio {
		input:
			cohort = cohort,
			container_registry = container_registry
	}

	Sample child = cohort.samples[parse_trio.child_index]
	Sample father = cohort.samples[parse_trio.father_index]
	Sample mother = cohort.samples[parse_trio.mother_index]

	scatter (sample in [child, father, mother]) {
		scatter (movie_bam in sample.movie_bams) {
			call SamtoolsFasta.samtools_fasta {
				input:
					bam = movie_bam.data,
					container_registry = container_registry
			}
		}
	}

	call yak_count as yak_count_father {
		input:
			sample_id = father.sample_id,
			reads_fastas = samtools_fasta.reads_fasta[1],
			container_registry = container_registry
	}

	call yak_count as yak_count_mother {
		input:
			sample_id = mother.sample_id,
			reads_fastas = samtools_fasta.reads_fasta[2],
			container_registry = container_registry
	}

	# Father is haplotype 1; mother is haplotype 2
	Map[String, String] haplotype_key_map = {
		"hap1": father.sample_id,
		"hap2": mother.sample_id
	}

	call AssembleGenome.assemble_genome {
		input:
			sample_id = cohort.cohort_id,
			reads_fastas = flatten(samtools_fasta.reads_fasta),
			reference = reference,
			hifiasm_extra_params = "-c1 -d1",
			father_yak = yak_count_father.yak,
			mother_yak = yak_count_mother.yak,
			container_registry = container_registry
	}

	scatter (zipped_assembly_fasta in assemble_genome.zipped_assembly_fastas) {
		call yak_trioeval {
			input:
				zipped_assembly_fasta = zipped_assembly_fasta,
				father_yak = yak_count_father.yak,
				mother_yak = yak_count_mother.yak,
				container_registry = container_registry
		}
	}

	output {
		Map[String, String] haplotype_key = haplotype_key_map
		Array[File] assembly_noseq_gfas = assemble_genome.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = assemble_genome.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = assemble_genome.zipped_assembly_fastas
		Array[File] assembly_stats = assemble_genome.assembly_stats
		IndexData asm_bam = assemble_genome.asm_bam
		Array[File] trioeval = yak_trioeval.trioeval
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		reference: {help: "Reference genome data"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task parse_trio {
	input {
		Cohort cohort

		String container_registry
	}

	command <<<
		set -euo pipefail

		parse_cohort.py \
			--cohort_json ~{write_json(cohort)} \
			--parse_trio
	>>>

	output {
		Int child_index = read_int("child_index.txt")
		Int father_index = read_int("father_index.txt")
		Int mother_index = read_int("mother_index.txt")
	}

	runtime {
		docker: "~{container_registry}/parse_cohort:0.0.1"
		cpu: 4
		memory: "14 GB"
		disk: "20 GB"
		preemptible: true
		maxRetries: 3
	}
}

task yak_count {
	input {
		String sample_id
		Array[File] reads_fastas

		String container_registry
	}

	Int threads = 32
	Int disk_size = ceil(size(reads_fastas[0], "GB") * length(reads_fastas) * 2 + 20)

	command <<<
		set -euo pipefail

		yak count \
			-t ~{threads} \
			-o ~{sample_id}.yak \
			~{sep=' ' reads_fastas}
	>>>

	output {
		File yak = "~{sample_id}.yak"
	}

	runtime {
		docker: "~{container_registry}/yak:b1a46c6"
		cpu: threads
		memory: "96 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task yak_trioeval {
	input {
		File zipped_assembly_fasta

		File father_yak
		File mother_yak

		String container_registry
	}

	String zipped_assembly_fasta_basename = basename(zipped_assembly_fasta, ".gz")
	Int threads = 16
	Int disk_size = ceil(size(zipped_assembly_fasta, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		yak trioeval \
			-t ~{threads} \
			~{father_yak} \
			~{mother_yak} \
			~{zipped_assembly_fasta} \
		> ~{zipped_assembly_fasta_basename}.trioeval.txt
	>>>

	output {
		File trioeval = "~{zipped_assembly_fasta_basename}.trioeval.txt"
	}

	runtime {
		docker: "~{container_registry}/yak:b1a46c6"
		cpu: threads
		memory: "48 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}
