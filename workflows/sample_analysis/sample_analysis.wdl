version 1.0

import "../common/structs.wdl"
import "../smrtcell_analysis/smrtcell_analysis.wdl" as SmrtcellAnalysis
import "../common/tasks/mosdepth.wdl" as Mosdepth
import "../common/tasks/pbsv_call.wdl" as PbsvCall
import "../common/tasks/zip_index_vcf.wdl" as ZipIndexVcf
import "../phase_vcf/phase_vcf.wdl" as PhaseVcf

workflow sample_analysis {
	input {
		Sample sample

		ReferenceData reference

		String deepvariant_version
		DeepVariantModel? deepvariant_model

		String container_registry
	}

	call SmrtcellAnalysis.smrtcell_analysis {
		input:
			sample = sample,
			reference = reference,
			container_registry = container_registry
	}

	scatter (bam_object in smrtcell_analysis.aligned_bams) {
		File aligned_bam = bam_object.data
		File aligned_bam_index = bam_object.data_index
	}

	# Run deepvariant_make_examples in parallel
	## The total number of tasks to split deepvariant_make_examples into
	Int deepvariant_total_num_tasks = 64
	## The number of workers to split tasks across; each worker will run (total_num_tasks / num_shards) tasks
	Int deepvariant_num_shards = 8
	Int deepvariant_num_tasks_per_shard = deepvariant_total_num_tasks / deepvariant_num_shards

	scatter (shard_index in range(deepvariant_num_shards)) {
		Int task_start_index = shard_index * deepvariant_num_tasks_per_shard
		Int task_end_index = (shard_index + 1) * deepvariant_num_tasks_per_shard - 1

		call deepvariant_make_examples {
			input:
				sample_id = sample.sample_id,
				aligned_bams = aligned_bam,
				aligned_bam_indices = aligned_bam_index,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				task_start_index = task_start_index,
				task_end_index = task_end_index,
				total_num_tasks = deepvariant_total_num_tasks,
				num_tasks_per_shard = deepvariant_num_tasks_per_shard,
				deepvariant_version = deepvariant_version
		}
	}

	call deepvariant_call_variants {
		input:
			sample_id = sample.sample_id,
			reference_name = reference.name,
			example_tfrecords = flatten(deepvariant_make_examples.example_tfrecords),
			deepvariant_model = deepvariant_model,
			deepvariant_version = deepvariant_version
	}

	call deepvariant_postprocess_variants {
		input:
			sample_id = sample.sample_id,
			tfrecord = deepvariant_call_variants.tfrecord,
			nonvariant_site_tfrecords = flatten(deepvariant_make_examples.nonvariant_site_tfrecords),
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			deepvariant_version = deepvariant_version
	}

	call bcftools_stats {
		input:
			vcf = deepvariant_postprocess_variants.vcf,
			params = "--apply-filters PASS --samples ~{sample.sample_id}",
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	call bcftools_roh {
		input:
			vcf = deepvariant_postprocess_variants.vcf,
			container_registry = container_registry
	}

	call PbsvCall.pbsv_call {
		input:
			sample_id = sample.sample_id,
			svsigs = smrtcell_analysis.svsigs,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			container_registry = container_registry
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = pbsv_call.pbsv_vcf,
			container_registry = container_registry
	}

	call PhaseVcf.phase_vcf {
		input:
			vcf = {"data": deepvariant_postprocess_variants.vcf, "data_index": deepvariant_postprocess_variants.vcf_index},
			aligned_bams = smrtcell_analysis.aligned_bams,
			reference = reference,
			container_registry = container_registry
	}

	scatter (bam_object in smrtcell_analysis.aligned_bams) {
		call whatshap_haplotag {
			input:
				phased_vcf = phase_vcf.phased_vcf.data,
				phased_vcf_index = phase_vcf.phased_vcf.data_index,
				aligned_bam = bam_object.data,
				aligned_bam_index = bam_object.data_index,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				container_registry = container_registry
		}
	}

	call merge_bams {
		input:
			bams = whatshap_haplotag.haplotagged_bam,
			output_bam_name = "~{sample.sample_id}.~{reference.name}.haplotagged.bam",
			container_registry = container_registry
	}

	call Mosdepth.mosdepth {
		input:
			aligned_bam = merge_bams.merged_bam,
			aligned_bam_index = merge_bams.merged_bam_index,
			container_registry = container_registry
	}

	call trgt {
		input:
			bam = merge_bams.merged_bam,
			bam_index = merge_bams.merged_bam_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			tandem_repeat_bed = reference.trgt_tandem_repeat_bed,
			container_registry = container_registry
	}

	call trgt_coverage_dropouts {
		input:
			bam = merge_bams.merged_bam,
			bam_index = merge_bams.merged_bam_index,
			output_prefix = "~{sample.sample_id}.~{reference.name}",
			tandem_repeat_bed = reference.trgt_tandem_repeat_bed,
			container_registry = container_registry
	}

	call cpg_pileup {
		input:
			bam = merge_bams.merged_bam,
			bam_index = merge_bams.merged_bam_index,
			output_prefix = "~{sample.sample_id}.~{reference.name}",
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			container_registry = container_registry
	}

	output {
		# smrtcell_analysis output
		Array[File] bam_stats = smrtcell_analysis.bam_stats
		Array[File] read_length_summary = smrtcell_analysis.read_length_summary
		Array[File] read_quality_summary = smrtcell_analysis.read_quality_summary
		Array[IndexData] aligned_bams = smrtcell_analysis.aligned_bams
		Array[File] aligned_bam_mosdepth_summary = smrtcell_analysis.aligned_bam_mosdepth_summary
		Array[File] aligned_bam_mosdepth_region_bed = smrtcell_analysis.aligned_bam_mosdepth_region_bed
		Array[File] svsigs = smrtcell_analysis.svsigs

		IndexData small_variant_vcf = {"data": deepvariant_postprocess_variants.vcf, "data_index": deepvariant_postprocess_variants.vcf_index}
		IndexData small_variant_gvcf = {"data": deepvariant_postprocess_variants.gvcf, "data_index": deepvariant_postprocess_variants.gvcf_index}
		File small_variant_vcf_stats = bcftools_stats.stats
		File small_variant_roh_bed = bcftools_roh.roh_bed
		IndexData sv_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		IndexData phased_small_variant_vcf = phase_vcf.phased_vcf
		File whatshap_stats_gtf = phase_vcf.whatshap_stats_gtf
		File whatshap_stats_tsv = phase_vcf.whatshap_stats_tsv
		File whatshap_stats_blocklist = phase_vcf.whatshap_stats_blocklist
		IndexData merged_haplotagged_bam = {"data": merge_bams.merged_bam, "data_index": merge_bams.merged_bam_index}
		File haplotagged_bam_mosdepth_summary = mosdepth.summary
		File haplotagged_bam_mosdepth_region_bed = mosdepth.region_bed
		IndexData trgt_spanning_reads = {"data": trgt.spanning_reads, "data_index": trgt.spanning_reads_index}
		IndexData trgt_repeat_vcf = {"data": trgt.repeat_vcf, "data_index": trgt.repeat_vcf_index}
		File trgt_dropouts = trgt_coverage_dropouts.trgt_dropouts
		Array[File] cpg_pileups = cpg_pileup.pileups
	}

	parameter_meta {
		sample: {help: "Sample information and associated data files"}
		reference: {help: "Reference genome data"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		container_registry: {help: "Container registry where docker images are hosted"}
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
		Int task_end_index
		Int total_num_tasks
		Int num_tasks_per_shard
		String deepvariant_version
	}

	Int disk_size = ceil(size(aligned_bams[0], "GB") * length(aligned_bams) * total_num_tasks * 2 + 400)

	command <<<
		set -euo pipefail

		seq ~{task_start_index} ~{task_end_index} \
		| parallel \
			--jobs ~{num_tasks_per_shard} \
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
				--examples ~{sample_id}.examples.tfrecord@~{total_num_tasks}.gz \
				--gvcf ~{sample_id}.gvcf.tfrecord@~{total_num_tasks}.gz \
				--task {}
	>>>

	output {
		Array[File] example_tfrecords = glob("~{sample_id}.examples.tfrecord*.gz")
		Array[File] nonvariant_site_tfrecords = glob("~{sample_id}.gvcf.tfrecord*.gz")
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: num_tasks_per_shard
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
		DeepVariantModel? deepvariant_model

		String deepvariant_version
	}

	String deepvariant_model_path = if (defined(deepvariant_model)) then sub(select_first([deepvariant_model]).model.data, "\\.data.*", "") else "/opt/models/pacbio/model.ckpt"
	Int disk_size = ceil(size(example_tfrecords[0], "GB") * length(example_tfrecords) * 2 + 200)

	command <<<
		set -euo pipefail

		/opt/deepvariant/bin/call_variants \
			--outfile ~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz \
			--examples ~{sep=' ' example_tfrecords} \
			--checkpoint ~{deepvariant_model_path}
	>>>

	output {
		File tfrecord = "~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: 64
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

		String deepvariant_version
	}

	Int disk_size = ceil((size(tfrecord, "GB") + size(reference, "GB") + size(nonvariant_site_tfrecords[0], "GB") * length(nonvariant_site_tfrecords)) * 2 + 20)

	command <<<
		set -euo pipefail

		/opt/deepvariant/bin/postprocess_variants \
			--ref ~{reference} \
			--infile ~{tfrecord} \
			--outfile ~{sample_id}.~{reference_name}.deepvariant.vcf.gz \
			--nonvariant_site_tfrecord_path ~{sep=' ' nonvariant_site_tfrecords} \
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
		cpu: 2
		memory: "32 GB"
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

task whatshap_haplotag {
	input {
		File phased_vcf
		File phased_vcf_index

		File aligned_bam
		File aligned_bam_index

		File reference
		File reference_index

		String container_registry
	}

	String bam_basename = basename(aligned_bam, ".bam")
	Int threads = 8
	Int disk_size = ceil((size(phased_vcf, "GB") + size(aligned_bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap haplotag \
			--tag-supplementary \
			--output-threads ~{threads} \
			--reference ~{reference} \
			--output ~{bam_basename}.haplotagged.bam \
			~{phased_vcf} \
			~{aligned_bam}
	>>>

	output {
		File haplotagged_bam = "~{bam_basename}.haplotagged.bam"
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
		File merged_bam_index = "~{output_bam_name}.bai"
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

task trgt {
	input {
		File bam
		File bam_index

		File reference
		File reference_index
		File tandem_repeat_bed

		String container_registry
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 4
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		trgt \
			--genome ~{reference} \
			--repeats ~{tandem_repeat_bed} \
			--reads ~{bam} \
			--output-prefix ~{bam_basename}.trgt

		bcftools sort \
			--output-type z \
			--output ~{bam_basename}.trgt.sorted.vcf.gz \
			~{bam_basename}.trgt.vcf.gz

		bcftools index \
			--tbi \
			~{bam_basename}.trgt.sorted.vcf.gz

		samtools sort \
			-@ ~{threads - 1} \
			-o ~{bam_basename}.trgt.spanning.sorted.bam \
			~{bam_basename}.trgt.spanning.bam

		samtools index \
			-@ ~{threads - 1} \
			~{bam_basename}.trgt.spanning.sorted.bam
	>>>

	output {
		File spanning_reads = "~{bam_basename}.trgt.spanning.sorted.bam"
		File spanning_reads_index = "~{bam_basename}.trgt.spanning.sorted.bam.bai"
		File repeat_vcf = "~{bam_basename}.trgt.sorted.vcf.gz"
		File repeat_vcf_index = "~{bam_basename}.trgt.sorted.vcf.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/trgt:v0.3.4"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task trgt_coverage_dropouts {
	input {
		File bam
		File bam_index

		String output_prefix

		File tandem_repeat_bed

		String container_registry
	}

	Int disk_size = ceil(size(bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		check_trgt_coverage.py \
			~{tandem_repeat_bed} \
			~{bam} \
		> ~{output_prefix}.trgt.dropouts.txt
	>>>

	output {
		File trgt_dropouts = "~{output_prefix}.trgt.dropouts.txt"
	}

	runtime {
		docker: "~{container_registry}/tandem-genotypes:07f9162"
		cpu: 4
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task cpg_pileup {
	input {
		File bam
		File bam_index

		String output_prefix

		File reference
		File reference_index

		String container_registry
	}

	Int threads = 48
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		/opt/scripts/pb-CpG-tools/aligned_bam_to_cpg_scores.py \
			--bam ~{bam} \
			--fasta ~{reference} \
			--output_label ~{output_prefix} \
			--threads ~{threads} \
			--min_mapq 1 \
			--modsites denovo \
			--pileup_mode model \
			--model_dir /opt/scripts/pb-CpG-tools/pileup_calling_model \
			--min_coverage 10
	>>>

	output {
		Array[File] pileups = glob("~{output_prefix}.{combined,hap1,hap2}.denovo.{bed,bw,mincov10.bed,mincov10.bw}")
	}

	runtime {
		docker: "~{container_registry}/pb-cpg-tools:b1a46c6"
		cpu: threads
		memory: "192 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}
