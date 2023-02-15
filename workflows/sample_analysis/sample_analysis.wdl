version 1.0

import "../common/structs.wdl"
import "../smrtcell_analysis/smrtcell_analysis.wdl" as SmrtcellAnalysis
import "../deepvariant/deepvariant.wdl" as DeepVariant
import "../common/tasks/bcftools_stats.wdl" as BcftoolsStats
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
		Boolean preemptible
	}

	call SmrtcellAnalysis.smrtcell_analysis {
		input:
			sample = sample,
			reference = reference,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call DeepVariant.deepvariant {
		input:
			sample_id = sample.sample_id,
			aligned_bams = smrtcell_analysis.aligned_bams,
			reference = reference,
			deepvariant_version = deepvariant_version,
			deepvariant_model = deepvariant_model,
			preemptible = preemptible
	}

	call BcftoolsStats.bcftools_stats {
		input:
			vcf = deepvariant.vcf.data,
			params = "--apply-filters PASS --samples ~{sample.sample_id}",
			reference = reference.fasta.data,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call bcftools_roh {
		input:
			vcf = deepvariant.vcf.data,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call PbsvCall.pbsv_call {
		input:
			sample_id = sample.sample_id,
			svsigs = smrtcell_analysis.svsigs,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = pbsv_call.pbsv_vcf,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call PhaseVcf.phase_vcf {
		input:
			vcf = deepvariant.vcf,
			aligned_bams = smrtcell_analysis.aligned_bams,
			reference = reference,
			container_registry = container_registry,
			preemptible = preemptible
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
				container_registry = container_registry,
				preemptible = preemptible
		}
	}

	call merge_bams {
		input:
			bams = whatshap_haplotag.haplotagged_bam,
			output_bam_name = "~{sample.sample_id}.~{reference.name}.haplotagged.bam",
			container_registry = container_registry,
			preemptible = preemptible
	}

	call Mosdepth.mosdepth {
		input:
			aligned_bam = merge_bams.merged_bam,
			aligned_bam_index = merge_bams.merged_bam_index,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call trgt {
		input:
			bam = merge_bams.merged_bam,
			bam_index = merge_bams.merged_bam_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			tandem_repeat_bed = reference.trgt_tandem_repeat_bed,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call trgt_coverage_dropouts {
		input:
			bam = merge_bams.merged_bam,
			bam_index = merge_bams.merged_bam_index,
			output_prefix = "~{sample.sample_id}.~{reference.name}",
			tandem_repeat_bed = reference.trgt_tandem_repeat_bed,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call cpg_pileup {
		input:
			bam = merge_bams.merged_bam,
			bam_index = merge_bams.merged_bam_index,
			output_prefix = "~{sample.sample_id}.~{reference.name}",
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			container_registry = container_registry,
			preemptible = preemptible
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

		IndexData small_variant_vcf = deepvariant.vcf
		IndexData small_variant_gvcf = deepvariant.gvcf
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

task bcftools_roh {
	input {
		File vcf

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/bcftools:1.14"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
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
		Boolean preemptible
	}

	String bam_basename = basename(aligned_bam, ".bam")
	Int threads = 4
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
		docker: "~{container_registry}/whatshap:1.4"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task merge_bams {
	input {
		Array[File] bams

		String output_bam_name

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/samtools:1.14"
		cpu: threads
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
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
		Boolean preemptible
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
		docker: "~{container_registry}/trgt:0.3.4"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
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
		Boolean preemptible
	}

	Int disk_size = ceil(size(bam, "GB") * 4 + 20)

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
		docker: "~{container_registry}/tandem-genotypes:1.8.1"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
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
		Boolean preemptible
	}

	Int threads = 48
	# Uses ~7 GB memory / thread
	Int mem_gb = threads * 8
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		aligned_bam_to_cpg_scores.py \
			--bam ~{bam} \
			--fasta ~{reference} \
			--output_label ~{output_prefix} \
			--threads ~{threads} \
			--min_mapq 1 \
			--modsites denovo \
			--pileup_mode model \
			--model_dir "$PILEUP_MODEL_DIR" \
			--min_coverage 10
	>>>

	output {
		Array[File] pileups = glob("~{output_prefix}.{combined,hap1,hap2}.denovo.{bed,bw,mincov10.bed,mincov10.bw}")
	}

	runtime {
		docker: "~{container_registry}/pb-cpg-tools:fed1a7b"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
