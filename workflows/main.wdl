version 1.0

import "common/structs.wdl"
import "smrtcell_analysis/smrtcell_analysis.wdl" as SmrtcellAnalysis
import "sample_analysis/sample_analysis.wdl" as SampleAnalysis
import "de_novo_assembly/de_novo_assembly.wdl" as DeNovoAssembly

workflow humanwgs {
	input {
		Sample sample

		IndexData reference_genome
		File reference_tandem_repeat_bed
		Array[String] chromosomes
		File reference_chromosome_lengths
		File trgt_tandem_repeat_bed

		String deepvariant_version
		File? deepvariant_model

		Boolean run_de_novo_assembly = false

		String container_registry
	}

	String reference_name = basename(reference_genome.data, ".fasta")

	call SmrtcellAnalysis.smrtcell_analysis {
		input:
			sample = sample,
			reference_genome = reference_genome,
			reference_name = reference_name,
			container_registry = container_registry
	}

	call SampleAnalysis.sample_analysis {
		input:
			sample_id = sample.sample_id,
			aligned_bams = smrtcell_analysis.aligned_bams,
			reference_genome = reference_genome,
			reference_name = reference_name,
			reference_tandem_repeat_bed = reference_tandem_repeat_bed,
			chromosomes = chromosomes,
			reference_chromosome_lengths = reference_chromosome_lengths,
			trgt_tandem_repeat_bed = trgt_tandem_repeat_bed,
			deepvariant_version = deepvariant_version,
			deepvariant_model = deepvariant_model,
			container_registry = container_registry
	}

	if (run_de_novo_assembly) {
		call DeNovoAssembly.de_novo_assembly {
			input:
				sample = sample,
				reference_genome = reference_genome,
				reference_name = reference_name,
				container_registry = container_registry
		}
	}

	output {
		# smrtcells_analysis output
		Array[File] bam_stats = smrtcell_analysis.bam_stats
		Array[File] read_length_summary = smrtcell_analysis.read_length_summary
		Array[File] read_quality_summary = smrtcell_analysis.read_quality_summary
		Array[IndexData] aligned_bams = smrtcell_analysis.aligned_bams
		Array[File] aligned_bam_mosdepth_global = smrtcell_analysis.aligned_bam_mosdepth_global
		Array[File] aligned_bam_mosdepth_region = smrtcell_analysis.aligned_bam_mosdepth_region
		Array[File] aligned_bam_mosdepth_summary = smrtcell_analysis.aligned_bam_mosdepth_summary
		Array[File] aligned_bam_mosdepth_region_bed = smrtcell_analysis.aligned_bam_mosdepth_region_bed

		# sample_analysis output
		IndexData pbsv_vcf = sample_analysis.pbsv_vcf
		IndexData deepvariant_vcf = sample_analysis.deepvariant_vcf
		IndexData deepvariant_gvcf = sample_analysis.deepvariant_gvcf
		File deepvariant_vcf_stats = sample_analysis.deepvariant_vcf_stats
		File deepvariant_roh_bed = sample_analysis.deepvariant_roh_bed
		IndexData phased_vcf = sample_analysis.phased_vcf
		File whatshap_stats_gtf = sample_analysis.whatshap_stats_gtf
		File whatshap_stats_tsv = sample_analysis.whatshap_stats_tsv
		File whatshap_stats_blocklist = sample_analysis.whatshap_stats_blocklist
		IndexData merged_haplotagged_bam = sample_analysis.merged_haplotagged_bam
		File haplotagged_bam_mosdepth_global = sample_analysis.haplotagged_bam_mosdepth_global
		File haplotagged_bam_mosdepth_region = sample_analysis.haplotagged_bam_mosdepth_region
		File haplotagged_bam_mosdepth_summary = sample_analysis.haplotagged_bam_mosdepth_summary
		File haplotagged_bam_mosdepth_region_bed = sample_analysis.haplotagged_bam_mosdepth_region_bed
		IndexData trgt_spanning_reads = sample_analysis.trgt_spanning_reads
		IndexData trgt_repeat_vcf = sample_analysis.trgt_repeat_vcf
		File trgt_dropouts = sample_analysis.trgt_dropouts
		Array[File] cpg_pileups = sample_analysis.cpg_pileups

		# de_novo_assembly output
		Array[File]? assembly_noseq_gfas = de_novo_assembly.assembly_noseq_gfas
		Array[File]? assembly_lowQ_beds = de_novo_assembly.assembly_lowQ_beds
		Array[File]? zipped_assembly_fastas = de_novo_assembly.zipped_assembly_fastas
		Array[File]? assembly_stats = de_novo_assembly.assembly_stats
		IndexData? asm_bam = de_novo_assembly.asm_bam
		IndexData? htsbox_vcf = de_novo_assembly.htsbox_vcf
		File? htsbox_vcf_stats = de_novo_assembly.htsbox_vcf_stats
	}

	parameter_meta {
		sample: {help: "Sample ID and unaligned movie bams and indices associated with the sample"}
		reference_genome: {help: "Reference genome and index to align reads to"}
		reference_tandem_repeat_bed: {help: "Tandem repeat locations in the reference genome"}
		chromosomes: {help: "Chromosomes to phase during WhatsHap phasing"}
		reference_chromosome_lengths: {help: "File specifying the lengths of each of the reference chromosomes"}
		trgt_tandem_repeat_bed: {help: "Repeat bed used by TRGT to output spanning reads and a repeat VCF"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		run_de_novo_assembly: {help: "Run the de novo assembly pipeline [false]"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}
