version 1.0

# Call variants using DeepVariant

import "../../structs.wdl"

workflow deepvariant {
	input {
		String sample_id
		Array[IndexData] aligned_bams

		IndexData reference_fasta
		String reference_name

		String deepvariant_version
		DeepVariantModel? deepvariant_model

		RuntimeAttributes default_runtime_attributes
	}

	scatter (bam_object in aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.data_index
	}

	Int total_deepvariant_tasks = 64
	Int num_shards = 8
	Int tasks_per_shard = total_deepvariant_tasks / num_shards

	scatter (shard_index in range(num_shards)) {
		Int task_start_index = shard_index * tasks_per_shard

		call deepvariant_make_examples {
			input:
				sample_id = sample_id,
				aligned_bams = aligned_bam,
				aligned_bam_indices = aligned_bam_index,
				reference = reference_fasta.data,
				reference_index = reference_fasta.data_index,
				task_start_index = task_start_index,
				tasks_per_shard = tasks_per_shard,
				total_deepvariant_tasks = total_deepvariant_tasks,
				deepvariant_version = deepvariant_version,
				runtime_attributes = default_runtime_attributes
		}
	}

	call deepvariant_call_variants {
		input:
			sample_id = sample_id,
			reference_name = reference_name,
			example_tfrecord_tars = deepvariant_make_examples.example_tfrecord_tar,
			deepvariant_model = deepvariant_model,
			total_deepvariant_tasks = total_deepvariant_tasks,
			deepvariant_version = deepvariant_version,
			runtime_attributes = default_runtime_attributes
	}

	call deepvariant_postprocess_variants {
		input:
			sample_id = sample_id,
			tfrecord = deepvariant_call_variants.tfrecord,
			nonvariant_site_tfrecord_tars = deepvariant_make_examples.nonvariant_site_tfrecord_tar,
			reference = reference_fasta.data,
			reference_index = reference_fasta.data_index,
			reference_name = reference_name,
			total_deepvariant_tasks = total_deepvariant_tasks,
			deepvariant_version = deepvariant_version,
			runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData vcf = {"data": deepvariant_postprocess_variants.vcf, "data_index": deepvariant_postprocess_variants.vcf_index}
		IndexData gvcf = {"data": deepvariant_postprocess_variants.gvcf, "data_index": deepvariant_postprocess_variants.gvcf_index}
	}

	parameter_meta {
		sample_id: {help: "Sample ID; used for naming files"}
		aligned_bams: {help: "Bam and index aligned to the reference genome for each movie associated with all samples in the cohort"}
		reference: {help: "Reference genome data"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task deepvariant_make_examples {
	input {
		String sample_id
		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference
		File reference_index

		Int task_start_index
		Int tasks_per_shard

		Int total_deepvariant_tasks
		String deepvariant_version

		RuntimeAttributes runtime_attributes
	}

	Int task_end_index = task_start_index + tasks_per_shard - 1
	Int disk_size = ceil(size(aligned_bams[0], "GB") * length(aligned_bams) * 2 + 50)
	Int mem_gb = tasks_per_shard * 4

	command <<<
		set -euo pipefail

		mkdir example_tfrecords nonvariant_site_tfrecords

		echo "DeepVariant version: $VERSION"

		seq ~{task_start_index} ~{task_end_index} \
		| parallel \
			--jobs ~{tasks_per_shard} \
			--halt 2 \
			/opt/deepvariant/bin/make_examples \
				--norealign_reads \
				--vsc_min_fraction_indels 0.12 \
				--pileup_image_width 199 \
				--track_ref_reads \
				--phase_reads \
				--partition_size=25000 \
				--max_reads_per_partition=600 \
				--alt_aligned_pileup=diff_channels \
				--add_hp_channel \
				--sort_by_haplotypes \
				--parse_sam_aux_fields \
				--min_mapping_quality=1 \
				--mode calling \
				--ref ~{reference} \
				--reads ~{sep="," aligned_bams} \
				--examples example_tfrecords/~{sample_id}.examples.tfrecord@~{total_deepvariant_tasks}.gz \
				--gvcf nonvariant_site_tfrecords/~{sample_id}.gvcf.tfrecord@~{total_deepvariant_tasks}.gz \
				--task {}

		tar -zcvf ~{sample_id}.~{task_start_index}.example_tfrecords.tar.gz example_tfrecords
		tar -zcvf ~{sample_id}.~{task_start_index}.nonvariant_site_tfrecords.tar.gz nonvariant_site_tfrecords
	>>>

	output {
		File example_tfrecord_tar = "~{sample_id}.~{task_start_index}.example_tfrecords.tar.gz"
		File nonvariant_site_tfrecord_tar = "~{sample_id}.~{task_start_index}.nonvariant_site_tfrecords.tar.gz"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: tasks_per_shard
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

task deepvariant_call_variants {
	input {
		String sample_id
		String reference_name
		Array[File] example_tfrecord_tars

		DeepVariantModel? deepvariant_model
		Int total_deepvariant_tasks
		String deepvariant_version

		RuntimeAttributes runtime_attributes
	}

	Int mem_gb = total_deepvariant_tasks * 4
	Int disk_size = ceil(size(example_tfrecord_tars[0], "GB") * length(example_tfrecord_tars) * 2 + 100)

	command <<<
		set -euo pipefail

		while read -r tfrecord_tar || [[ -n "${tfrecord_tar}" ]]; do
			tar -zxvf "${tfrecord_tar}"
		done < ~{write_lines(example_tfrecord_tars)}

		deepvariant_model_path=~{if (defined(deepvariant_model)) then sub(select_first([deepvariant_model]).model.data, "\\.data.*", "") else "/opt/models/pacbio/model.ckpt"}

		echo "DeepVariant version: $VERSION"

		/opt/deepvariant/bin/call_variants \
			--outfile ~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz \
			--examples "example_tfrecords/~{sample_id}.examples.tfrecord@~{total_deepvariant_tasks}.gz" \
			--checkpoint "${deepvariant_model_path}"
	>>>

	output {
		File tfrecord = "~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: total_deepvariant_tasks
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

task deepvariant_postprocess_variants {
	input {
		String sample_id
		File tfrecord
		Array[File] nonvariant_site_tfrecord_tars

		File reference
		File reference_index
		String reference_name

		Int total_deepvariant_tasks
		String deepvariant_version

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(tfrecord, "GB") + size(reference, "GB") + size(nonvariant_site_tfrecord_tars[0], "GB") * length(nonvariant_site_tfrecord_tars)) * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r nonvariant_site_tfrecord_tar || [[ -n "${nonvariant_site_tfrecord_tar}" ]]; do
			tar -zxvf "${nonvariant_site_tfrecord_tar}"
		done < ~{write_lines(nonvariant_site_tfrecord_tars)}

		echo "DeepVariant version: $VERSION"

		/opt/deepvariant/bin/postprocess_variants \
			--vcf_stats_report=false \
			--ref ~{reference} \
			--infile ~{tfrecord} \
			--outfile ~{sample_id}.~{reference_name}.deepvariant.vcf.gz \
			--nonvariant_site_tfrecord_path "nonvariant_site_tfrecords/~{sample_id}.gvcf.tfrecord@~{total_deepvariant_tasks}.gz" \
			--gvcf_outfile ~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz
	>>>

	output {
		File vcf = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz"
		File vcf_index = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz.tbi"
		File gvcf = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz"
		File gvcf_index = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz.tbi"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: 2
		memory: "32 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
