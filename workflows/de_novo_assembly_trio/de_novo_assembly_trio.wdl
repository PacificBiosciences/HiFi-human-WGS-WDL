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

	call parse_families {
		input:
			cohort_json = write_json(cohort),
			container_registry = container_registry
	}

	# Run de_novo_assembly for each child with mother and father samples present in the cohort
	# Multiple children per family and multiple unrelated families can be included in the cohort and will each produce separate child assemblies
	scatter (family in parse_families.families) {
		Sample father = cohort.samples[family.father_index]
		Sample mother = cohort.samples[family.mother_index]

		scatter (movie_bam in father.movie_bams) {
			call SamtoolsFasta.samtools_fasta as samtools_fasta_father {
				input:
					bam = movie_bam,
					container_registry = container_registry
			}
		}

		call yak_count as yak_count_father {
			input:
				sample_id = father.sample_id,
				reads_fastas = samtools_fasta_father.reads_fasta,
				container_registry = container_registry
		}

		scatter (movie_bam in mother.movie_bams) {
			call SamtoolsFasta.samtools_fasta as samtools_fasta_mother {
				input:
					bam = movie_bam,
					container_registry = container_registry
			}
		}

		call yak_count as yak_count_mother {
			input:
				sample_id = mother.sample_id,
				reads_fastas = samtools_fasta_mother.reads_fasta,
				container_registry = container_registry
		}

		# Father is haplotype 1; mother is haplotype 2
		Map[String, String] haplotype_key_map = {
			"hap1": father.sample_id,
			"hap2": mother.sample_id
		}

		scatter (child_index in family.child_indices) {
			Sample child = cohort.samples[child_index]

			scatter (movie_bam in child.movie_bams) {
				call SamtoolsFasta.samtools_fasta as samtools_fasta_child {
					input:
						bam = movie_bam,
						container_registry = container_registry
				}
			}

			call AssembleGenome.assemble_genome {
				input:
					sample_id = "~{cohort.cohort_id}.~{child.sample_id}",
					reads_fastas = flatten([
							samtools_fasta_child.reads_fasta,
							samtools_fasta_father.reads_fasta,
							samtools_fasta_mother.reads_fasta
						]),
					reference = reference,
					hifiasm_extra_params = "-c1 -d1",
					father_yak = yak_count_father.yak,
					mother_yak = yak_count_mother.yak,
					container_registry = container_registry
			}
		}
	}

	output {
		Array[Map[String, String]] haplotype_key = haplotype_key_map
		Array[Array[File]] assembly_noseq_gfas = flatten(assemble_genome.assembly_noseq_gfas)
		Array[Array[File]] assembly_lowQ_beds = flatten(assemble_genome.assembly_lowQ_beds)
		Array[Array[File]] zipped_assembly_fastas = flatten(assemble_genome.zipped_assembly_fastas)
		Array[Array[File]] assembly_stats = flatten(assemble_genome.assembly_stats)
		Array[IndexData] asm_bams = flatten(assemble_genome.asm_bam)
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		reference: {help: "Reference genome data"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task parse_families {
	input {
		File cohort_json

		String container_registry
	}

	command <<<
		set -euo pipefail

		parse_cohort.py \
			--cohort_json ~{cohort_json} \
			--parse_families
	>>>

	output {
		Array[FamilySampleIndices] families = read_json("families.json")
	}

	runtime {
		docker: "~{container_registry}/parse_cohort:1.0.0"
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
