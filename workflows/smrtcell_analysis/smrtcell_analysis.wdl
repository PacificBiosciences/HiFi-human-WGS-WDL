version 1.0

struct Sample {
	String sample_id
	Array[File] movie_bams
}

workflow smrtcell_analysis {
	input {
		Sample sample

		File reference
		File reference_index

		String container_registry
	}

	scatter (unaligned_bam in sample.movie_bams) {
		call smrtcell_stats {
			input:
				sample_id = sample.sample_id,
				bam = unaligned_bam,
				container_registry = container_registry
		}

		call pbmm2_align {
			input:
				sample_id = sample.sample_id,
				bam = unaligned_bam,
				reference = reference,
				reference_index = reference_index,
				container_registry = container_registry
		}

		call mosdepth {
			input:
				aligned_bam = pbmm2_align.aligned_bam,
				aligned_bam_index = pbmm2_align.aligned_bam_index,
				container_registry = container_registry
		}
	}

	output {
		Array[File] bam_stats = smrtcell_stats.bam_stats
		Array[File] read_length_summary = smrtcell_stats.read_length_summary
		Array[File] read_quality_summary = smrtcell_stats.read_quality_summary
		Array[File] aligned_bam = pbmm2_align.aligned_bam
		Array[File] aligned_bam_index = pbmm2_align.aligned_bam_index
		Array[File] mosdepth_global = mosdepth.global
		Array[File] mosdepth_region = mosdepth.region
		Array[File] mosdepth_summary = mosdepth.summary
		Array[File] mosdepth_region_bed = mosdepth.region_bed
	}

	parameter_meta {
		sample: {help: "Sample ID"}
		unaligned_bam: {help: "Unaligned movie BAM output by the sequencer"}
		reference: {help: "Reference genome to align reads to"}
		reference_index: {help: "Index for the reference genome"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task smrtcell_stats {
	input {
		String sample_id
		File bam

		String container_registry
	}

	String movie = basename(bam, ".bam")
	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		extract_read_length_and_qual.py \
			~{bam} \
		> ~{sample_id}.~{movie}.read_length_and_quality.tsv

		awk '{{ b=int($2/1000); b=(b>39?39:b); print 1000*b "\t" $2; }}' \
			~{sample_id}.~{movie}.read_length_and_quality.tsv \
			| sort -k1,1g \
			| datamash -g 1 count 1 sum 2 \
			| awk 'BEGIN {{ for(i=0;i<=39;i++) {{ print 1000*i"\t0\t0"; }} }} {{ print; }}' \
			| sort -k1,1g \
			| datamash -g 1 sum 2 sum 3 \
		> ~{sample_id}.~{movie}.read_length_summary.tsv

		awk '{{ print ($3>50?50:$3) "\t" $2; }}' \
				~{sample_id}.~{movie}.read_length_and_quality.tsv \
			| sort -k1,1g \
			| datamash -g 1 count 1 sum 2 \
			| awk 'BEGIN {{ for(i=0;i<=60;i++) {{ print i"\t0\t0"; }} }} {{ print; }}' \
			| sort -k1,1g \
			| datamash -g 1 sum 2 sum 3 \
		> ~{sample_id}.~{movie}.read_quality_summary.tsv
	>>>

	output {
		File bam_stats = "~{sample_id}.~{movie}.read_length_and_quality.tsv"
		File read_length_summary = "~{sample_id}.~{movie}.read_length_summary.tsv"
		File read_quality_summary = "~{sample_id}.~{movie}.read_quality_summary.tsv"
	}

	runtime {
		docker: "~{container_registry}/smrtcell_stats:b1a46c6"
		cpu: 4
		memory: "24 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task pbmm2_align {
	input {
		String sample_id
		File bam

		File reference
		File reference_index

		String container_registry
	}

	String movie = basename(bam, ".bam")
	String ref_name = basename(reference, ".fasta")

	Int threads = 24
	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		pbmm2 align \
			--num-threads ~{threads} \
			--preset CCS \
			--sample ~{sample_id} \
			--log-level INFO \
			--sort \
			--unmapped \
			-c 0 \
			-y 70 \
			~{reference} \
			~{bam} \
			~{sample_id}.~{movie}.~{ref_name}.aligned.bam
	>>>

	output {
		File aligned_bam = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
		File aligned_bam_index = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
	}

	runtime {
		docker: "~{container_registry}/pbmm2:b1a46c6"
		cpu: threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task mosdepth {
	input {
		File aligned_bam
		File aligned_bam_index

		String container_registry
	}

	String prefix = basename(aligned_bam, ".bam")
	Int threads = 4
	Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

	command <<<
		mosdepth \
			--threads ~{threads} \
			--by 500 \
			--no-per-base \
			--use-median \
			~{prefix} \
			~{aligned_bam}
	>>>

	output {
		File global = "~{prefix}.mosdepth.global.dist.txt"
		File region = "~{prefix}.mosdepth.region.dist.txt"
		File summary = "~{prefix}.mosdepth.summary.txt"
		File region_bed = "~{prefix}.regions.bed.gz"
	}

	runtime {
		docker: "~{container_registry}/mosdepth:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}
