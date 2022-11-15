version 1.0

import "../common/structs.wdl"
import "../common/common.wdl" as common

workflow sample_analysis {
	input {
		String sample_id
		Array[IndexData] aligned_bams

		IndexData reference_genome
		String reference_name
		File reference_tandem_repeat_bed
		Array[String] chromosomes
		File reference_chromosome_lengths

		String deepvariant_version
		File? deepvariant_model
		String container_registry
	}

	Int deepvariant_threads = 64

	scatter (bam_object in aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.data_index

		call pbsv_discover {
			input:
				aligned_bam = aligned_bam,
				aligned_bam_index = aligned_bam_index,
				reference = reference_genome.data,
				reference_index = reference_genome.data_index,
				reference_tandem_repeat_bed = reference_tandem_repeat_bed,
				container_registry = container_registry
		}
	}

	call pbsv_call {
		input:
			sample_id = sample_id,
			svsigs = pbsv_discover.svsig,
			reference = reference_genome.data,
			reference_name = reference_name,
			container_registry = container_registry
	}

	call deepvariant_make_examples {
		input:
			sample_id = sample_id,
			aligned_bams = aligned_bam,
			aligned_bam_indices = aligned_bam_index,
			reference = reference_genome.data,
			reference_index = reference_genome.data_index,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version
	}

	call deepvariant_call_variants {
		input:
			sample_id = sample_id,
			reference_name = reference_name,
			example_tfrecords = deepvariant_make_examples.example_tfrecords,
			deepvariant_model = deepvariant_model,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version
	}

	call deepvariant_postprocess_variants {
		input:
			sample_id = sample_id,
			tfrecord = deepvariant_call_variants.tfrecord,
			nonvariant_site_tfrecords = deepvariant_make_examples.nonvariant_site_tfrecords,
			reference = reference_genome.data,
			reference_index = reference_genome.data_index,
			reference_name = reference_name,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version
	}

	call bcftools_stats {
		input:
			vcf = deepvariant_postprocess_variants.vcf,
			params = "--apply-filters PASS --samples ~{sample_id}",
			reference = reference_genome.data,
			container_registry = container_registry
	}

	call bcftools_roh {
		input:
			vcf = deepvariant_postprocess_variants.vcf,
			container_registry = container_registry
	}

	scatter (chromosome in chromosomes) {
		call split_vcf {
			input:
				vcf = deepvariant_postprocess_variants.vcf,
				region = chromosome,
				container_registry = container_registry
		}

		call whatshap_phase {
			input:
				vcf = split_vcf.region_vcf,
				vcf_index = split_vcf.region_vcf_index,
				chromosome = chromosome,
				aligned_bams = aligned_bam,
				aligned_bam_indices = aligned_bam_index,
				reference = reference_genome.data,
				container_registry = container_registry
		}
	}

	call bcftools_concat {
		input:
			vcfs = whatshap_phase.phased_vcf,
			vcf_indices = whatshap_phase.phased_vcf_index,
			output_vcf_name = "~{sample_id}.~{reference_name}.deepvariant.phased.vcf.gz",
			container_registry = container_registry
	}

	call whatshap_stats {
		input:
			phased_vcf = bcftools_concat.concatenated_vcf,
			phased_vcf_index = bcftools_concat.concatenated_vcf_index,
			reference_chromosome_lengths = reference_chromosome_lengths,
			container_registry = container_registry
	}

	scatter (bam_object in aligned_bams) {
		call whatshap_haplotag {
			input:
				phased_vcf = bcftools_concat.concatenated_vcf,
				phased_vcf_index = bcftools_concat.concatenated_vcf_index,
				aligned_bam = bam_object.data,
				aligned_bam_index = bam_object.data_index,
				output_bam_name = basename(bam_object.data, ".aligned.bam") + ".deepvariant.haplotagged.bam",
				reference = reference_genome.data,
				container_registry = container_registry
		}
	}

	call merge_bams {
		input:
			bams = whatshap_haplotag.haplotagged_bam,
			output_bam_name = "~{sample_id}.~{reference_name}.deepvariant.haplotagged.bam",
			container_registry = container_registry
	}

	call common.mosdepth {
		input:
			aligned_bam = merge_bams.merged_bam,
			aligned_bam_index = merge_bams.merged_bam_index,
			container_registry = container_registry
	}

	output {
		File pbsv_vcf = pbsv_call.pbsv_vcf
		IndexData deepvariant_vcf = {"data": deepvariant_postprocess_variants.vcf, "data_index": deepvariant_postprocess_variants.vcf_index}
		IndexData deepvariant_gvcf = {"data": deepvariant_postprocess_variants.gvcf, "data_index": deepvariant_postprocess_variants.gvcf_index}
		File deepvariant_vcf_stats = bcftools_stats.stats
		File deepvariant_roh_bed = bcftools_roh.roh_bed
		IndexData phased_vcf = {"data": bcftools_concat.concatenated_vcf, "data_index": bcftools_concat.concatenated_vcf_index}
		File whatshap_stats_gtf = whatshap_stats.gtf
		File whatshap_stats_tsv = whatshap_stats.tsv
		File whatshap_stats_blocklist = whatshap_stats.blocklist
		IndexData merged_haplotagged_bam = {"data": merge_bams.merged_bam, "data_index": merge_bams.merged_bam_index}
		File haplotagged_bam_mosdepth_global = mosdepth.global
		File haplotagged_bam_mosdepth_region = mosdepth.region
		File haplotagged_bam_mosdepth_summary = mosdepth.summary
		File haplotagged_bam_mosdepth_region_bed = mosdepth.region_bed
	}

	parameter_meta {
		sample_id: {help: "Sample ID"}
		aligned_bams: {help: "Bam and index aligned to the reference genome for each movie associated with the sample"}
		reference_genome: {help: "Reference genome and index to align reads to"}
		reference_name: {help: "Basename of the reference genome; used for file naming"}
		reference_tandem_repeat_bed: {help: "Tandem repeat locations in the reference genome"}
		chromosomes: {help: "Chromosomes to phase during WhatsHap phasing"}
		reference_chromosome_lengths: {help: "File specifying the lengths of each of the reference chromosomes"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task pbsv_discover {
	input {
		File aligned_bam
		File aligned_bam_index

		File reference
		File reference_index
		File reference_tandem_repeat_bed

		String container_registry
	}

	String prefix = basename(aligned_bam, ".bam")
	Int disk_size = ceil((size(aligned_bam, "GB") + size(reference, "GB") + size(reference_tandem_repeat_bed, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv discover \
			--log-level INFO \
			--hifi \
			--tandem-repeats ~{reference_tandem_repeat_bed} \
			~{aligned_bam} \
			~{prefix}.svsig.gz
	>>>

	output {
		File svsig = "~{prefix}.svsig.gz"
	}

	runtime {
		docker: "~{container_registry}/pbsv:b1a46c6"
		cpu: 4
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task pbsv_call {
	input {
		String sample_id
		Array[File] svsigs

		File reference
		String reference_name

		String container_registry
	}

	Int threads = 8
	Int disk_size = ceil((size(svsigs[0], "GB") * length(svsigs) + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv call \
			--hifi \
			-m 20 \
			--num-threads ~{threads} \
			~{reference} \
			~{sep=' ' svsigs} \
			~{sample_id}.~{reference_name}.pbsv.vcf
	>>>

	output {
		File pbsv_vcf = "~{sample_id}.~{reference_name}.pbsv.vcf"
	}

	runtime {
		docker: "~{container_registry}/pbsv:b1a46c6"
		cpu: threads
		memory: "48 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
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
	}

	Int disk_size = 500

	command <<<
		set -euo pipefail

		seq 0 ~{deepvariant_threads} \
		| parallel \
			--jobs ~{deepvariant_threads} \
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
		Array[File] example_tfrecords = glob("~{sample_id}.examples.tfrecord.*.gz")
		Array[File] nonvariant_site_tfrecords = glob("~{sample_id}.gvcf.tfrecord.*.gz")
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: deepvariant_threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task deepvariant_call_variants {
	input {
		String sample_id
		String reference_name
		Array[File] example_tfrecords
		File? deepvariant_model

		Int deepvariant_threads
		String deepvariant_version
	}

	String example_tfrecord_path = sub(example_tfrecords[0], "/" + basename(example_tfrecords[0]), "")
	Int disk_size = 500

	command <<<
		set -euo pipefail

		if [[ -n "~{deepvariant_model}" ]]; then
			DEEPVARIANT_MODEL="~{deepvariant_model}"
		else
			DEEPVARIANT_MODEL=/opt/models/pacbio/model.ckpt
		fi

		/opt/deepvariant/bin/call_variants \
			--outfile ~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz \
			--examples ~{example_tfrecord_path}/~{sample_id}.examples.tfrecord@~{deepvariant_threads}.gz \
			--checkpoint "$DEEPVARIANT_MODEL"
	>>>

	output {
		File tfrecord = "~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: deepvariant_threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
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
	}

	String nonvariant_site_tfrecord_path = sub(nonvariant_site_tfrecords[0], "/" + basename(nonvariant_site_tfrecords[0]), "")
	Int disk_size = ceil((size(tfrecord, "GB") + size(reference, "GB") + size(nonvariant_site_tfrecords[0], "GB") * length(nonvariant_site_tfrecords)) * 2 + 20)

	command <<<
		set -euo pipefail

		/opt/deepvariant/bin/postprocess_variants \
			--ref ~{reference} \
			--infile ~{tfrecord} \
			--outfile ~{sample_id}.~{reference_name}.deepvariant.vcf.gz \
			--nonvariant_site_tfrecord_path ~{nonvariant_site_tfrecord_path}/~{sample_id}.gvcf.tfrecord@~{deepvariant_threads}.gz \
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
		cpu: 4
		memory: "30 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcftools_stats {
	input {
		File vcf
		String? params

		File? reference

		String container_registry
	}

	String vcf_basename = basename(vcf, ".gz")

	Int threads = 4
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools stats \
			--threads ~{threads - 1} \
			~{params} \
			~{"--fasta-ref " + reference} \
			~{vcf} \
		> ~{vcf_basename}.stats.txt
	>>>

	output {
		File stats = "~{vcf_basename}.stats.txt"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcftools_roh {
	input {
		File vcf

		String container_registry
	}

	String vcf_basename = basename(vcf, ".vcf.gz")

	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		echo -e "#chr\tstart\tend\tqual" > ~{vcf_basename}.roh.bed
		bcftools roh \
			--AF-dflt 0.4 \
			~{vcf} \
		| awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' \
		>> ~{vcf_basename}.roh.bed
	>>>

	output {
		File roh_bed = "~{vcf_basename}.roh.bed"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task split_vcf {
	input {
		File vcf
		String region

		String container_registry
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		tabix \
			-h \
			~{vcf} \
			~{region} \
		> ~{vcf_basename}.~{region}.vcf

		bgzip ~{vcf_basename}.~{region}.vcf
		tabix ~{vcf_basename}.~{region}.vcf.gz
	>>>

	output {
		File region_vcf = "~{vcf_basename}.~{region}.vcf.gz"
		File region_vcf_index = "~{vcf_basename}.~{region}.vcf.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/htslib:b1a46c6"
		cpu: 4
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task whatshap_phase {
	input {
		File vcf
		File vcf_index
		String chromosome

		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference

		String container_registry
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(aligned_bams[0], "GB") * length(aligned_bams)) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap phase \
			--indels \
			--reference ~{reference} \
			--chromosome ~{chromosome} \
			--output ~{vcf_basename}.phased.vcf.gz \
			~{vcf} \
			~{sep=' ' aligned_bams}

		tabix ~{vcf_basename}.phased.vcf.gz
	>>>

	output {
		File phased_vcf = "~{vcf_basename}.phased.vcf.gz"
		File phased_vcf_index = "~{vcf_basename}.phased.vcf.gz.tbi"
	}

	# TODO doesn't whatshap phase run single-threaded? why the giant machine?
	runtime {
		docker: "~{container_registry}/whatshap:b1a46c6"
		cpu: 32
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcftools_concat {
	input {
		Array[File] vcfs
		Array[File] vcf_indices
		String output_vcf_name

		String container_registry
	}

	Int disk_size = ceil(size(vcfs[0], "GB") * length(vcfs) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools concat \
			--allow-overlaps \
			--output ~{output_vcf_name} \
			--output-type z \
			~{sep=' ' vcfs}

		tabix "~{output_vcf_name}"
	>>>

	output {
		File concatenated_vcf = "~{output_vcf_name}"
		File concatenated_vcf_index = "~{output_vcf_name}.tbi"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: 4
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task whatshap_stats {
	input {
		File phased_vcf
		File phased_vcf_index

		File reference_chromosome_lengths

		String container_registry
	}

	String output_basename = basename(phased_vcf, ".vcf.gz")
	Int disk_size = ceil(size(phased_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap stats \
			--gtf ~{output_basename}.gtf \
			--tsv ~{output_basename}.tsv \
			--block-list ~{output_basename}.blocklist \
			--chr-lengths ~{reference_chromosome_lengths} \
			~{phased_vcf}
	>>>

	output {
		File gtf = "~{output_basename}.gtf"
		File tsv = "~{output_basename}.tsv"
		File blocklist = "~{output_basename}.blocklist"
	}

	runtime {
		docker: "~{container_registry}/whatshap:b1a46c6"
		cpu: 4
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task whatshap_haplotag {
	input {
		File phased_vcf
		File phased_vcf_index

		File aligned_bam
		File aligned_bam_index

		String output_bam_name

		File reference

		String container_registry
	}

	Int threads = 8
	Int disk_size = ceil((size(phased_vcf, "GB") + size(aligned_bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap haplotag \
			--tag-supplementary \
			--output-threads ~{threads} \
			--reference ~{reference} \
			--output ~{output_bam_name} \
			~{phased_vcf} \
			~{aligned_bam}
	>>>

	output {
		File haplotagged_bam = "~{output_bam_name}"
	}

	runtime {
		docker: "~{container_registry}/whatshap:b1a46c6"
		cpu: threads
		memory: "30 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task merge_bams {
	input {
		Array[File] bams

		String output_bam_name

		String container_registry
	}

	Int threads = 8
	Int disk_size = ceil(size(bams[0], "GB") * length(bams) * 2 + 20)

	command <<<
		set -euo pipefail

		if [[ "~{length(bams)}" -eq 1 ]]; then
			mv ~{bams[0]} ~{output_bam_name}
		else
			samtools merge \
				-@ ~{threads - 1} \
				-o ~{output_bam_name} \
				~{sep=' ' bams}
		fi

		samtools index ~{output_bam_name}
	>>>

	output {
		File merged_bam = "~{output_bam_name}"
		File merged_bam_index = "~{output_bam_name}.tbi"
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
