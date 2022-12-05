version 1.0

import "../common/structs.wdl"
import "../common/common.wdl" as common
import "../common/phase_vcf.wdl" as PhaseVCF

workflow cohort_analysis {
	input {
		String cohort_id
		File pedigree
		Array[IndexData] aligned_bams
		Array[File] svsigs
		Array[IndexData] gvcfs

		ReferenceData reference

		File slivar_js

		String container_registry
	}

	scatter (gvcf_object in gvcfs) {
		File gvcf = gvcf_object.data
		File gvcf_index = gvcf_object.data_index
	}

	call common.pbsv_call {
		input:
			sample_id = cohort_id,
			svsigs = svsigs,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			container_registry = container_registry
	}

	call common.zip_index_vcf {
		input:
			vcf = pbsv_call.pbsv_vcf,
			container_registry = container_registry
	}

	call glnexus {
		input:
			cohort_id = cohort_id,
			gvcfs = gvcf,
			gvcf_indices = gvcf_index,
			reference_name = reference.name
	}

	call bcf_to_vcf {
		input:
			bcf = glnexus.bcf,
			container_registry = container_registry
	}

	call PhaseVCF.phase_vcf {
		input:
			vcf = {"data": bcf_to_vcf.vcf, "data_index": bcf_to_vcf.vcf_index},
			aligned_bams = aligned_bams,
			reference = reference,
			container_registry = container_registry
	}

	call bcftools_norm {
		input:
			vcf = phase_vcf.phased_vcf.data,
			vcf_index = phase_vcf.phased_vcf.data_index,
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	call slivar_small_variant {
		input:
			bcf = bcftools_norm.normalized_bcf,
			bcf_index = bcftools_norm.normalized_bcf_index,
			pedigree = pedigree,
			reference = reference.fasta.data,
			slivar_js = slivar_js,
			gnomad_af = reference.gnomad_af,
			hprc_af = reference.hprc_af,
			gff = reference.gff,
			container_registry = container_registry
	}

	output {
		IndexData sv_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		IndexData phased_joint_called_vcf = phase_vcf.phased_vcf
		File whatshap_stats_gtf = phase_vcf.whatshap_stats_gtf
		File whatshap_stats_tsv = phase_vcf.whatshap_stats_tsv
		File whatshap_stats_blocklist = phase_vcf.whatshap_stats_blocklist
	}

	parameter_meta {
		cohort_id: {help: "Sample ID for the cohort; used for naming files"}
		pedigree: {help: "Pedigree for the cohort"}
		aligned_bams: {help: "Bam and index aligned to the reference genome for each movie associated with the sample"}
		svsigs: {help: "pbsv svsig files for each sample and movie bam in the cohort"}
		gvcfs: {help: "gVCFs and indices for each sample in the cohort"}
		reference: {help: "ReferenceData"}
		slivar_js: {help: "Additional javascript functions for slivar"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task glnexus {
	input {
		String cohort_id
		Array[File] gvcfs
		Array[File] gvcf_indices

		String reference_name
	}

	Int threads = 24
	Int mem_gbytes = 30
	Int disk_size = ceil((size(gvcfs[0], "GB") * length(gvcfs)) * 2 + 20)

	command <<<
		glnexus_cli \
			--threads ~{threads} \
			--mem-gbytes ~{mem_gbytes} \
			--dir ~{cohort_id}.~{reference_name}.GLnexus.DB \
			--config DeepVariant_unfiltered \
			~{sep=' ' gvcfs} \
		> ~{cohort_id}.~{reference_name}.deepvariant.glnexus.bcf
	>>>

	output {
		File bcf = "~{cohort_id}.~{reference_name}.deepvariant.glnexus.bcf"
	}

	runtime {
		docker: "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
		cpu: threads
		memory: mem_gbytes + " GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcf_to_vcf {
	input {
		File bcf

		String container_registry
	}

	String bcf_basename = basename(bcf, ".bcf")
	Int threads = 4
	Int disk_size = ceil(size(bcf, "GB") * 2 + 20)

	command <<<
		bcftools view \
			--threads ~{threads} \
			--output-type z \
			--output ~{bcf_basename}.vcf.gz \
			~{bcf}

		bcftoosl index ~{bcf_basename}.vcf.gz
	>>>

	output {
		File vcf = "~{bcf_basename}.vcf.gz"
		File vcf_index = "~{bcf_basename}.vcf.gz.tbi"
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

task bcftools_norm {
	input {
		File vcf
		File vcf_index

		File reference

		String container_registry
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int threads = 4
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		bcftools norm \
			--multiallelics \
			- \
			--output-type b \
			--fasta-ref ~{reference} \
			~{vcf} \
		| bcftools sort \
			--threads ~{threads} \
			--output-type b \
			-o ~{vcf_basename}.norm.bcf

		bcftools index ~{vcf_basename}.norm.bcf
	>>>

	output {
		File normalized_bcf = "~{vcf_basename}.norm.bcf"
		File normalized_bcf_index = "~{vcf_basename}.norm.bcf.csi"
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

task slivar_small_variant {
	input {
		File bcf
		File bcf_index
		File pedigree

		File reference

		File slivar_js
		File gnomad_af
		File hprc_af
		File gff

		String container_registry
	}

	String bcf_basename = basename(bcf, ".bcf")
	Int threads = 12

	command <<<
		pslivar \
			--process ~{threads} \
			--pass-only \
			--js ~{slivar_js} \
			# TODO filters
			--gnotate ~{gnomad_af} \
			--gnotate ~{hprc_af} \
			--vcf ~{bcf} \
			--ped ~{pedigree} \
		| bcftools csq \
			--local-csq \
			--samples - \
			--ncsq 40 \
			--gff-annot ~{gff} \
			--fasta-ref ~{reference} \
			- \
			--output-type z \
			--output ~{bcf_basename}.slivar.vcf.gz

		tabix ~{bcf_basename}.slivar.vcf.gz
	>>>

	output {
		File slivar_vcf = "~{bcf_basename}.slivar.vcf.gz"
		File slivar_vcf_index = "~{bcf_basename}.slivar.vcf.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/slivar:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: "500 GB"
		preemptible: true
		maxRetries: 3
	}
}
