version 1.0

# Assemble a genome using hifiasm. Can be used for single-sample or trio-based assembly.

import "../humanwgs_structs.wdl"

workflow assemble_genome {
	input {
		String sample_id
		Array[File] reads_fastas

		ReferenceData reference

		String? hifiasm_extra_params
		File? father_yak
		File? mother_yak

		String backend
		RuntimeAttributes default_runtime_attributes
		RuntimeAttributes on_demand_runtime_attributes
	}

	# Preemptible jobs in GCP cannot live >24h; force on-demand for assembly in GCP
	RuntimeAttributes assemble_runtime_attributes = if (backend == "GCP") then on_demand_runtime_attributes else default_runtime_attributes

	call hifiasm_assemble {
		input:
			sample_id = sample_id,
			reads_fastas = reads_fastas,
			extra_params = hifiasm_extra_params,
			father_yak = father_yak,
			mother_yak = mother_yak,
			runtime_attributes = assemble_runtime_attributes
	}

	scatter (gfa in hifiasm_assemble.assembly_hap_gfas) {
		# Convert gfa to zipped fasta; calculate assembly stats
		call gfa2fa {
			input:
				gfa = gfa,
				reference_index = reference.fasta.data_index,
				runtime_attributes = default_runtime_attributes
		}
	}

	call align_hifiasm {
		input:
			sample_id = sample_id,
			query_sequences = gfa2fa.zipped_fasta,
			reference = reference.fasta.data,
			reference_name = reference.name,
			runtime_attributes = default_runtime_attributes
	}

	output {
		Array[File] assembly_noseq_gfas = hifiasm_assemble.assembly_noseq_gfas
		Array[File] assembly_lowQ_beds = hifiasm_assemble.assembly_lowQ_beds
		Array[File] zipped_assembly_fastas = gfa2fa.zipped_fasta
		Array[File] assembly_stats = gfa2fa.assembly_stats
		IndexData asm_bam = {"data": align_hifiasm.asm_bam, "data_index": align_hifiasm.asm_bam_index}
	}

	parameter_meta {
		sample_id: {help: "Sample ID; used for naming files"}
		reads_fastas: {help: "Reads in fasta format to be used for assembly; one for each movie bam to be used in assembly. Reads fastas from one or more sample may be combined to use in the assembly"}
		reference: {help: "Reference genome data"}
		hiiasm_extra_params: {help: "[OPTIONAL] Additional parameters to pass to hifiasm assembly"}
		father_yak: {help: "[OPTIONAL] kmer counts for the father; required if running trio-based assembly"}
		mother_yak: {help: "[OPTIONAL] kmer counts for the mother; required if running trio-based assembly"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
		on_demand_runtime_attributes: {help: "RuntimeAttributes for tasks that require dedicated instances"}
	}
}

# Note that this task will run ~25% faster on intel vs. AMD processors
task hifiasm_assemble {
	input {
		String sample_id
		Array[File] reads_fastas

		String? extra_params
		File? father_yak
		File? mother_yak

		RuntimeAttributes runtime_attributes
	}

	String prefix = "~{sample_id}.asm"
	Int threads = 48
	Int mem_gb = threads * 6
	Int disk_size = ceil((size(reads_fastas[0], "GB") * length(reads_fastas)) * 4 + 20)

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
		docker: "~{runtime_attributes.container_registry}/hifiasm:0.19.4"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task gfa2fa {
	input {
		File gfa

		File reference_index

		RuntimeAttributes runtime_attributes
	}

	String gfa_basename = basename(gfa, ".gfa")
	Int disk_size = ceil(size(gfa, "GB") * 3 + 20)
	Int threads = 2

	command <<<
		set -euo pipefail

		gfatools gfa2fa \
			~{gfa} \
		> ~{gfa_basename}.fasta

		bgzip \
			--threads ~{threads} \
			--stdout \
			~{gfa_basename}.fasta \
		> ~{gfa_basename}.fasta.gz

		# Calculate assembly stats
		k8 \
			/opt/calN50/calN50.js \
			-f ~{reference_index} \
			~{gfa_basename}.fasta.gz \
		> ~{gfa_basename}.fasta.stats.txt
	>>>

	output {
		File zipped_fasta = "~{gfa_basename}.fasta.gz"
		File assembly_stats = "~{gfa_basename}.fasta.stats.txt"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/gfatools:0.4_01091f2"
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

task align_hifiasm {
	input {
		String sample_id
		Array[File] query_sequences

		File reference
		String reference_name

		RuntimeAttributes runtime_attributes
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
			~{reference} \
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
		docker: "~{runtime_attributes.container_registry}/align_hifiasm:2.17_1.14"
		cpu: threads
		memory: "256 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
