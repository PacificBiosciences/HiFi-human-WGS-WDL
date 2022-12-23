version 1.0

import "../common/structs.wdl"

workflow deepvariant {
	input {
		String sample_id
		Array[IndexData] aligned_bams
		ReferenceData reference
		String deepvariant_version
		DeepVariantModel? deepvariant_model
		Boolean preemptible
	}

	Int deepvariant_threads = 64

	scatter (bam_object in aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.data_index
	}

	call deepvariant_make_examples {
		input:
			sample_id = sample_id,
			aligned_bams = aligned_bam,
			aligned_bam_indices = aligned_bam_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version,
			preemptible = preemptible
	}

	call deepvariant_call_variants {
		input:
			sample_id = sample_id,
			reference_name = reference.name,
			example_tfrecords = deepvariant_make_examples.example_tfrecords,
			deepvariant_model = deepvariant_model,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version,
			preemptible = preemptible
	}

	call deepvariant_postprocess_variants {
		input:
			sample_id = sample_id,
			tfrecord = deepvariant_call_variants.tfrecord,
			nonvariant_site_tfrecords = deepvariant_make_examples.nonvariant_site_tfrecords,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version,
			preemptible = preemptible
	}

	output {
		IndexData vcf = {"data": deepvariant_postprocess_variants.vcf, "data_index": deepvariant_postprocess_variants.vcf_index}
		IndexData gvcf = {"data": deepvariant_postprocess_variants.gvcf, "data_index": deepvariant_postprocess_variants.gvcf_index}
	}
}


task deepvariant_make_examples {
	input {
		String sample_id
		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference
		File reference_index

		Int deepvariant_threads
		String deepvariant_version
		Boolean preemptible
	}

	Int disk_size = ceil(size(aligned_bams[0], "GB") * length(aligned_bams) * 2 + 50)
	Int mem_gb = deepvariant_threads * 4

	command <<<
		set -euo pipefail

		seq 0 ~{deepvariant_threads - 1} \
		| parallel \
			--jobs ~{deepvariant_threads} \
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
				--examples ~{sample_id}.examples.tfrecord@~{deepvariant_threads}.gz \
				--gvcf ~{sample_id}.gvcf.tfrecord@~{deepvariant_threads}.gz \
				--task {}
	>>>

	output {
		Array[File] example_tfrecords = glob("~{sample_id}.examples.tfrecord*.gz")
		Array[File] nonvariant_site_tfrecords = glob("~{sample_id}.gvcf.tfrecord*.gz")
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: deepvariant_threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task deepvariant_call_variants {
	input {
		String sample_id
		String reference_name
		Array[File] example_tfrecords

		DeepVariantModel? deepvariant_model
		Int deepvariant_threads
		String deepvariant_version
		Boolean preemptible
	}

	String deepvariant_model_path = if (defined(deepvariant_model)) then sub(select_first([deepvariant_model]).model.data, "\\.data.*", "") else "/opt/models/pacbio/model.ckpt"

	Int disk_size = ceil(size(example_tfrecords[0], "GB") * length(example_tfrecords) * 2 + 100)
	Int mem_gb = deepvariant_threads * 4

	command <<<
		set -euo pipefail

		# extract the path where the first example_tfrecord is located; all example_tfrecords will be located at the same base path
		example_tfrecord_dir=$(dirname ~{example_tfrecords[0]})

		/opt/deepvariant/bin/call_variants \
			--outfile ~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz \
			--examples "$example_tfrecord_dir/~{sample_id}.examples.tfrecord@~{deepvariant_threads}.gz" \
			--checkpoint ~{deepvariant_model_path}
	>>>

	output {
		File tfrecord = "~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: deepvariant_threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task deepvariant_postprocess_variants {
	input {
		String sample_id
		File tfrecord
		Array[File] nonvariant_site_tfrecords

		File reference
		File reference_index
		String reference_name

		Int deepvariant_threads
		String deepvariant_version
		Boolean preemptible
	}

	Int disk_size = ceil((size(tfrecord, "GB") + size(reference, "GB") + size(nonvariant_site_tfrecords[0], "GB") * length(nonvariant_site_tfrecords)) * 2 + 20)

	command <<<
		set -euo pipefail

		# extract the path where the first nonvariant_site_tfrecord is located; all nonvariant_site_tfrecord will be located at the same base path
		nonvariant_site_tfrecord_dir=$(dirname ~{nonvariant_site_tfrecords[0]})

		/opt/deepvariant/bin/postprocess_variants \
			--ref ~{reference} \
			--infile ~{tfrecord} \
			--outfile ~{sample_id}.~{reference_name}.deepvariant.vcf.gz \
			--nonvariant_site_tfrecord_path "$nonvariant_site_tfrecord_dir/~{sample_id}.gvcf.tfrecord@~{deepvariant_threads}.gz" \
			--gvcf_outfile ~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz
	>>>

	output {
		File vcf = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz"
		File vcf_index = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz.tbi"
		File gvcf = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz"
		File gvcf_index = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz.tbi"
		File report = "~{sample_id}.~{reference_name}.deepvariant.visual_report.html"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: 1
		memory: "32 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
