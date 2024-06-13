version 1.0

import "../../structs.wdl"

workflow hiphase {
	input {
		String id
		String refname

		Array[String] sample_ids
		Array[IndexData] vcfs
		Array[IndexData] bams

		Boolean haplotag = false

		IndexData reference_fasta

		RuntimeAttributes default_runtime_attributes
	}

	scatter (vcf_object in vcfs) {
		File vcf = vcf_object.data
		File vcf_index = vcf_object.data_index
		# generate an array of phased VCF names that match the input VCFs
		String phased_vcf_name = basename(vcf, ".vcf.gz") + ".phased.vcf.gz"
		String phased_vcf_index_name = basename(vcf, ".vcf.gz") + ".phased.vcf.gz.tbi"
	}

	scatter (bam_object in bams) {
		File bam = bam_object.data
		File bam_index = bam_object.data_index
	}

	# for samples, haplotag BAMs
	if (haplotag) {
		scatter (bam_object in bams) {
			# generate an array of haplotagged BAM names that match the input BAMs
			String haplotagged_bam_name = basename(bam_object.data, ".bam") + ".haplotagged.bam"
			String haplotagged_bam_index_name = basename(bam_object.data, ".bam") + ".haplotagged.bam.bai"
		}
	}

	call run_hiphase {
		input:
			id = id,
			refname = refname,
			sample_ids = sample_ids,
			vcfs = vcf,
			vcf_indices = vcf_index,
			phased_vcf_names = phased_vcf_name,
			phased_vcf_index_names = phased_vcf_index_name,
			bams = bam,
			bam_indices = bam_index,
			haplotagged_bam_names = select_first([haplotagged_bam_name, []]),
			haplotagged_bam_index_names = select_first([haplotagged_bam_index_name, []]),
			reference = reference_fasta.data,
			reference_index = reference_fasta.data_index,
			runtime_attributes = default_runtime_attributes
	}

	# zip phased VCFs and indices
	scatter (pair in zip(run_hiphase.phased_vcfs, run_hiphase.phased_vcf_indices)) {
		IndexData phased_vcf = {"data": pair.left, "data_index": pair.right}
	}

	# zip haplotagged BAMs and indices
	scatter (pair in zip(run_hiphase.haplotagged_bams, run_hiphase.haplotagged_bam_indices)) {
		IndexData haplotagged_bam = {"data": pair.left, "data_index": pair.right}
	}

	output {
		Array[IndexData] phased_vcfs = phased_vcf
		Array[IndexData] haplotagged_bams = haplotagged_bam
		File hiphase_stats = run_hiphase.hiphase_stats
		File hiphase_blocks = run_hiphase.hiphase_blocks
		File? hiphase_haplotags = run_hiphase.hiphase_haplotags
	}

	parameter_meta {
		id: {help: "Sample ID or Cohort ID"}
		sample_ids: {help: "Samples to phase"}
		haplotag: {help: "Optional: Whether to haplotag BAMs"}
		vcfs: {help: "VCFs to phase"}
		bams: {help: "Bam and index aligned to the reference genome for each movie associated with the sample"}
		refname: {help: "Reference genome short name"}
		reference: {help: "Reference genome fasta and index"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task run_hiphase {
	input {
		String id
		String refname

		Array[String] sample_ids

		Array[File] vcfs
		Array[File] vcf_indices
		Array[String] phased_vcf_names
		Array[String] phased_vcf_index_names

		Array[File] bams
		Array[File] bam_indices

		Array[String] haplotagged_bam_names
		Array[String] haplotagged_bam_index_names # !UnusedDeclaration; used along with haplotagged_bam_names

		File reference
		File reference_index

		RuntimeAttributes runtime_attributes
	}

	# to handle single samples with very high depth, we had to increase memory to 3*gpu
	# to handle cohorts with very high depth, we had to increase memory to 4*gpu
	Int threads = 16
	Int mem_gb = threads * 4
	Int disk_size = ceil(size(vcfs, "GB") + size(reference, "GB") + size(bams, "GB") * 2 + 20)
	String haplotags_param = if length(haplotagged_bam_names) > 0 then "--haplotag-file ~{id}.~{refname}.hiphase.haplotags.tsv" else ""

	command <<<
		set -euo pipefail

		hiphase --version

		# phase VCFs and haplotag BAM
		hiphase --threads ~{threads} \
			~{sep=" " prefix("--sample-name ", sample_ids)} \
			~{sep=" " prefix("--vcf ", vcfs)} \
			~{sep=" " prefix("--output-vcf ", phased_vcf_names)} \
			~{sep=" " prefix("--bam ", bams)} \
			~{true="--output-bam " false="" length(haplotagged_bam_names) > 0} ~{sep=" --output-bam " haplotagged_bam_names} \
			--reference ~{reference} \
			--summary-file ~{id}.~{refname}.hiphase.stats.tsv \
			--blocks-file ~{id}.~{refname}.hiphase.blocks.tsv \
			~{haplotags_param} \
			--global-realignment-cputime 300
	>>>

	output {
		Array[File] phased_vcfs = phased_vcf_names
		Array[File] phased_vcf_indices = phased_vcf_index_names
		Array[File] haplotagged_bams = haplotagged_bam_names
		Array[File] haplotagged_bam_indices = haplotagged_bam_index_names
		File hiphase_stats = "~{id}.~{refname}.hiphase.stats.tsv"
		File hiphase_blocks = "~{id}.~{refname}.hiphase.blocks.tsv"
		File? hiphase_haplotags = "~{id}.~{refname}.hiphase.haplotags.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hiphase@sha256:493ed4244608f29d7e2e180af23b20879c71ae3692201a610c7f1f980ee094e8"
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
