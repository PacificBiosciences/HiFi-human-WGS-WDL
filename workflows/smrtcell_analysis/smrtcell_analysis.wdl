version 1.0

import "../common/structs.wdl"
import "../common/common.wdl" as common

workflow smrtcell_analysis {
	input {
		Sample sample

		IndexData reference_genome
		String reference_name

		String container_registry
	}

	scatter (movie_bam in sample.movie_bams) {
		call smrtcell_stats {
			input:
				sample_id = sample.sample_id,
				bam = movie_bam.data,
				container_registry = container_registry
		}

		call pbmm2_align {
			input:
				sample_id = sample.sample_id,
				bam = movie_bam.data,
				reference = reference_genome.data,
				reference_index = reference_genome.data_index,
				reference_name = reference_name,
				container_registry = container_registry
		}

		call common.mosdepth {
			input:
				aligned_bam = pbmm2_align.aligned_bam,
				aligned_bam_index = pbmm2_align.aligned_bam_index,
				container_registry = container_registry
		}

		IndexData aligned_bam = {
			"data": pbmm2_align.aligned_bam,
			"data_index": pbmm2_align.aligned_bam_index
		}
	}

	output {
		Array[File] bam_stats = smrtcell_stats.bam_stats
		Array[File] read_length_summary = smrtcell_stats.read_length_summary
		Array[File] read_quality_summary = smrtcell_stats.read_quality_summary
		Array[IndexData] aligned_bams = aligned_bam
		Array[File] aligned_bam_mosdepth_global = mosdepth.global
		Array[File] aligned_bam_mosdepth_region = mosdepth.region
		Array[File] aligned_bam_mosdepth_summary = mosdepth.summary
		Array[File] aligned_bam_mosdepth_region_bed = mosdepth.region_bed
	}

	parameter_meta {
		sample: {help: "Sample ID and unaligned movie bams and indices associated with the sample"}
		reference_genome: {help: "Reference genome and index to align reads to"}
		reference_name: {help: "Basename of the reference genome; used for file naming"}
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
		String reference_name

		String container_registry
	}

	String movie = basename(bam, ".bam")

	Int threads = 24
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

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
			~{sample_id}.~{movie}.~{reference_name}.aligned.bam
	>>>

	output {
		File aligned_bam = "~{sample_id}.~{movie}.~{reference_name}.aligned.bam"
		File aligned_bam_index = "~{sample_id}.~{movie}.~{reference_name}.aligned.bam.bai"
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
