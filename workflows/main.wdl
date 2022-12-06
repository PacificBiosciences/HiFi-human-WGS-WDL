version 1.0

import "common/structs.wdl"
import "smrtcell_analysis/smrtcell_analysis.wdl" as SmrtcellAnalysis
import "sample_analysis/sample_analysis.wdl" as SampleAnalysis
import "de_novo_assembly/de_novo_assembly.wdl" as DeNovoAssembly
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis
import "slivar/slivar.wdl" as Slivar

workflow humanwgs {
	input {
		Cohort cohort

		# TODO host these and import the default values so users don't need to define them?
		ReferenceData reference

		# TODO organize into annotations struct?
		File slivar_js
		HpoData hpo
		File ensembl_to_hgnc
		File clinvar_lookup
		File lof_lookup

		String deepvariant_version
		File? deepvariant_model

		Boolean run_de_novo_assembly = false

		String container_registry
	}

	scatter (sample in cohort.samples) {
		call SmrtcellAnalysis.smrtcell_analysis {
			input:
				sample = sample,
				reference = reference,
				deepvariant_version = deepvariant_version,
				deepvariant_model = deepvariant_model,
				container_registry = container_registry
		}

		# TODO if this is run for every sample, move it into smrtcellAnalysis?
		call SampleAnalysis.sample_analysis {
			input:
				sample_id = sample.sample_id,
				small_variant_vcf = smrtcell_analysis.deepvariant_vcf,
				aligned_bams = smrtcell_analysis.aligned_bams,
				svsigs = smrtcell_analysis.svsigs,
				reference = reference,
				container_registry = container_registry
		}
	}

	if (length(cohort.samples) == 1) {
		if (run_de_novo_assembly) {
			call DeNovoAssembly.de_novo_assembly {
				input:
					sample = cohort.samples[0],
					reference = reference,
					container_registry = container_registry
			}
		}
	}

	if (length(cohort.samples) > 1) {
		call CohortAnalysis.cohort_analysis {
			input:
				cohort_id = cohort.cohort_id,
				aligned_bams = flatten(smrtcell_analysis.aligned_bams),
				svsigs = flatten(smrtcell_analysis.svsigs),
				gvcfs = smrtcell_analysis.deepvariant_gvcf,
				reference = reference,
				container_registry = container_registry
		}
	}

	IndexData slivar_input_vcf = select_first([cohort_analysis.phased_joint_called_vcf, sample_analysis.phased_small_variant_vcf[0]])

	call Slivar.slivar {
		input:
			cohort = cohort,
			vcf = slivar_input_vcf,
			reference = reference,
			hpo = hpo,
			ensembl_to_hgnc = ensembl_to_hgnc,
			slivar_js = slivar_js,
			lof_lookup = lof_lookup,
			clinvar_lookup = clinvar_lookup,
			container_registry = container_registry
	}

	output {
		# smrtcells_analysis output
		Array[Array[File]] bam_stats = smrtcell_analysis.bam_stats
		Array[Array[File]] read_length_summary = smrtcell_analysis.read_length_summary
		Array[Array[File]] read_quality_summary = smrtcell_analysis.read_quality_summary
		Array[Array[IndexData]] aligned_bams = smrtcell_analysis.aligned_bams
		Array[Array[File]] aligned_bam_mosdepth_global = smrtcell_analysis.aligned_bam_mosdepth_global
		Array[Array[File]] aligned_bam_mosdepth_region = smrtcell_analysis.aligned_bam_mosdepth_region
		Array[Array[File]] aligned_bam_mosdepth_summary = smrtcell_analysis.aligned_bam_mosdepth_summary
		Array[Array[File]] aligned_bam_mosdepth_region_bed = smrtcell_analysis.aligned_bam_mosdepth_region_bed
		Array[IndexData] deepvariant_vcfs = smrtcell_analysis.deepvariant_vcf
		Array[IndexData] deepvariant_gvcf = smrtcell_analysis.deepvariant_gvcf
		Array[File] deepvariant_vcf_stats = smrtcell_analysis.deepvariant_vcf_stats
		Array[File] deepvariant_roh_bed = smrtcell_analysis.deepvariant_roh_bed

		# sample_analysis output
		Array[IndexData] merged_haplotagged_bam = sample_analysis.merged_haplotagged_bam
		Array[File] haplotagged_bam_mosdepth_global = sample_analysis.haplotagged_bam_mosdepth_global
		Array[File] haplotagged_bam_mosdepth_region = sample_analysis.haplotagged_bam_mosdepth_region
		Array[File] haplotagged_bam_mosdepth_summary = sample_analysis.haplotagged_bam_mosdepth_summary
		Array[File] haplotagged_bam_mosdepth_region_bed = sample_analysis.haplotagged_bam_mosdepth_region_bed
		Array[IndexData] trgt_spanning_reads = sample_analysis.trgt_spanning_reads
		Array[IndexData] trgt_repeat_vcf = sample_analysis.trgt_repeat_vcf
		Array[File] trgt_dropouts = sample_analysis.trgt_dropouts
		Array[Array[File]] cpg_pileups = sample_analysis.cpg_pileups

		# output by sample_analysis and cohort_analysis
		## TODO might separate these into sample_ and cohort_ outputs for clarity
		Array[IndexData] sv_vcfs = select_all(flatten([sample_analysis.sv_vcf, [cohort_analysis.sv_vcf]]))
		Array[IndexData] phased_small_variant_vcfs = select_all(flatten([sample_analysis.phased_small_variant_vcf, [cohort_analysis.phased_joint_called_vcf]]))
		Array[File] whatshap_stats_gtfs = select_all(flatten([sample_analysis.whatshap_stats_gtf, [cohort_analysis.whatshap_stats_gtf]]))
		Array[File] whatshap_stats_tsvs = select_all(flatten([sample_analysis.whatshap_stats_tsv, [cohort_analysis.whatshap_stats_tsv]]))
		Array[File] whatshap_stats_blocklists = select_all(flatten([sample_analysis.whatshap_stats_blocklist, [cohort_analysis.whatshap_stats_blocklist]]))

		# slivar output
		IndexData filterd_vcf = slivar.filterd_vcf
		IndexData compound_het_vcf = slivar.compound_het_vcf
		File filtered_tsv = slivar.filtered_tsv
		File compound_het_tsv = slivar.compound_het_tsv

		# singleton de_novo_assembly output
		Array[File]? assembly_noseq_gfas = de_novo_assembly.assembly_noseq_gfas
		Array[File]? assembly_lowQ_beds = de_novo_assembly.assembly_lowQ_beds
		Array[File]? zipped_assembly_fastas = de_novo_assembly.zipped_assembly_fastas
		Array[File]? assembly_stats = de_novo_assembly.assembly_stats
		IndexData? asm_bam = de_novo_assembly.asm_bam
		IndexData? htsbox_vcf = de_novo_assembly.htsbox_vcf
		File? htsbox_vcf_stats = de_novo_assembly.htsbox_vcf_stats
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		reference: {help: "Reference genome data"}
		slivar_js: {help: "Additional javascript functions for slivar"}
		hpo: {help: "HPO annotation lookups"}
		ensembl_to_hgnc: {help: "Ensembl to HGNC gene mapping"}
		lof_lookup: {help: "LOF lookup for slivar annotations"}
		clinvar_lookup: {help: "ClinVar lookup for slivar annotations"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		run_de_novo_assembly: {help: "Run the de novo assembly pipeline [false]"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}
