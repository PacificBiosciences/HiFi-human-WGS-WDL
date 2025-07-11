version 1.0

import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "upstream/upstream.wdl" as Upstream
import "joint/joint.wdl" as Joint
import "downstream/downstream.wdl" as Downstream
import "wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "tertiary/tertiary.wdl" as TertiaryAnalysis
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities

workflow humanwgs_family {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline, with joint calling for related samples."
  }

  parameter_meta {
    family: {
      name: "Family struct describing samples, relationships, and unaligned BAM paths"
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    ref_map_file: {
      name: "TSV containing reference genome file paths; must match backend"
    }
    tertiary_map_file: {
      name: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    max_reads_per_alignment_chunk: {
      name: "Maximum reads per alignment chunk"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    glnexus_mem_gb: {
      name: "Override GLnexus memory request (GB)"
    }
    gpu: {
      name: "Use GPU when possible"
    }
    backend: {
      name: "Backend where the workflow will be executed",
      choices: ["GCP", "Azure", "AWS-HealthOmics", "HPC"]
    }
    zones: {
      name: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      help: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      name: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      name: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'"
    }
    preemptible: {
      name: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      name: "Debug version for testing purposes"
    }
  }

  input {
    Family family

    String phenotypes = "HP:0000001"

    File ref_map_file
    File? tertiary_map_file

    Int max_reads_per_alignment_chunk = 500000
    Int pharmcat_min_coverage = 10
    Int? glnexus_mem_gb

    Boolean gpu = false

    # Backend configuration
    String backend
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry

    Boolean preemptible = true

    String? debug_version
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = backend,
      zones              = zones,
      cpuPlatform        = cpuPlatform,
      gpuType            = gpuType,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  Map [String, String] ref_map = read_map(ref_map_file)

  Boolean single_sample = length(family.samples) == 1

  Map[String, String] pedigree_sex = {
    "MALE": "1",
    "FEMALE": "2",
    "": "."
  }

  scatter (sample in family.samples) {
    String sample_id = sample.sample_id
    Boolean is_trio_kid = defined(sample.father_id) && defined(sample.mother_id)  # !UnusedDeclaration
    Boolean is_duo_kid = defined(sample.father_id) != defined(sample.mother_id)   # !UnusedDeclaration

    call Upstream.upstream {
      input:
        sample_id                     = sample.sample_id,
        sex                           = sample.sex,
        hifi_reads                    = sample.hifi_reads,
        ref_map_file                  = ref_map_file,
        max_reads_per_alignment_chunk = max_reads_per_alignment_chunk,
        single_sample                 = single_sample,
        gpu                           = gpu,
        default_runtime_attributes    = default_runtime_attributes
    }

    # write sample metadata similar to pedigree format
    # family_id, sample_id, father_id, mother_id, sex, affected
    Array[String] sample_metadata = [
      family.family_id,
      sample.sample_id,
      select_first([sample.father_id, "."]),
      select_first([sample.mother_id, "."]),
      pedigree_sex[upstream.inferred_sex],
      if sample.affected then "2" else "1"
    ]
  }

  if (!single_sample) {
    call Joint.joint {
      input:
        family_id                  = family.family_id,
        sample_ids                 = sample_id,
        gvcfs                      = upstream.small_variant_gvcf,
        gvcf_indices               = upstream.small_variant_gvcf_index,
        discover_tars              = upstream.discover_tar,
        aligned_bams               = upstream.out_bam,
        aligned_bam_indices        = upstream.out_bam_index,
        ref_map_file               = ref_map_file,
        glnexus_mem_gb             = glnexus_mem_gb,
        default_runtime_attributes = default_runtime_attributes
    }
  }

  scatter (sample_index in range(length(family.samples))) {
    call Downstream.downstream {
      input:
        sample_id                  = sample_id[sample_index],
        small_variant_vcf          = select_first([joint.split_joint_small_variant_vcfs, upstream.small_variant_vcf])[sample_index],
        small_variant_vcf_index    = select_first([joint.split_joint_small_variant_vcf_indices, upstream.small_variant_vcf_index])[sample_index],
        sv_vcf                     = select_first([joint.split_joint_structural_variant_vcfs, select_all(upstream.sv_vcf)])[sample_index],
        sv_vcf_index               = select_first([joint.split_joint_structural_variant_vcf_indices, select_all(upstream.sv_vcf_index)])[sample_index],
        trgt_vcf                   = upstream.trgt_vcf[sample_index],
        trgt_vcf_index             = upstream.trgt_vcf_index[sample_index],
        aligned_bam                = upstream.out_bam[sample_index],
        aligned_bam_index          = upstream.out_bam_index[sample_index],
        pharmcat_min_coverage      = pharmcat_min_coverage,
        ref_map_file               = ref_map_file,
        default_runtime_attributes = default_runtime_attributes
    }
  }

  if (!single_sample) {
    call Bcftools.bcftools_merge as merge_small_variant_vcfs {
      input:
        vcfs               = downstream.phased_small_variant_vcf,
        vcf_indices        = downstream.phased_small_variant_vcf_index,
        out_prefix         = "~{family.family_id}.joint.~{ref_map['name']}.small_variants.phased",
        runtime_attributes = default_runtime_attributes
    }

    call Bcftools.bcftools_merge as merge_sv_vcfs {
      input:
        vcfs               = downstream.phased_sv_vcf,
        vcf_indices        = downstream.phased_sv_vcf_index,
        out_prefix         = "~{family.family_id}.joint.~{ref_map['name']}.structural_variants.phased",
        runtime_attributes = default_runtime_attributes
    }

    call Trgt.trgt_merge {
      input:
        vcfs               = downstream.phased_trgt_vcf,
        vcf_indices        = downstream.phased_trgt_vcf_index,
        ref_fasta          = ref_map["fasta"],                              # !FileCoercion
        ref_index          = ref_map["fasta_index"],                        # !FileCoercion
        out_prefix         = "~{family.family_id}.merged.~{ref_map['name']}.trgt",
        runtime_attributes = default_runtime_attributes
    }
  }

  if (defined(tertiary_map_file)) {
    call TertiaryAnalysis.tertiary_analysis {
      input:
        sample_metadata            = sample_metadata,
        phenotypes                 = phenotypes,
        is_trio_kid                = is_trio_kid,
        is_duo_kid                 = is_duo_kid,
        small_variant_vcf          = select_first([merge_small_variant_vcfs.merged_vcf, downstream.phased_small_variant_vcf[0]]),
        small_variant_vcf_index    = select_first([merge_small_variant_vcfs.merged_vcf_index, downstream.phased_small_variant_vcf_index[0]]),
        sv_vcf                     = select_first([merge_sv_vcfs.merged_vcf, downstream.phased_sv_vcf[0]]),
        sv_vcf_index               = select_first([merge_sv_vcfs.merged_vcf_index, downstream.phased_sv_vcf_index[0]]),
        ref_map_file               = ref_map_file,
        tertiary_map_file          = select_first([tertiary_map_file]),
        default_runtime_attributes = default_runtime_attributes
    }
  }

    Map[String, Array[String]] stats = {
    'sample_id': sample_id,
    'num_reads': downstream.stat_num_reads,
    'read_length_mean': downstream.stat_read_length_mean,
    'read_length_median': downstream.stat_read_length_median,
    'read_quality_mean': downstream.stat_read_quality_mean,
    'read_quality_median': downstream.stat_read_quality_median,
    'mapped_read_count': downstream.stat_mapped_read_count,
    'mapped_percent': downstream.stat_mapped_percent,
    'mean_depth': upstream.stat_mean_depth,
    'inferred_sex': upstream.inferred_sex,
    'stat_phased_basepairs': downstream.stat_phased_basepairs,
    'phase_block_ng50': downstream.stat_phase_block_ng50,
    'cpg_combined_count': downstream.stat_combined_cpg_count,
    'cpg_hap1_count': downstream.stat_hap1_cpg_count,
    'cpg_hap2_count': downstream.stat_hap2_cpg_count,
    'SNV_count': downstream.stat_SNV_count,
    'TSTV_ratio': downstream.stat_TSTV_ratio,
    'HETHOM_ratio': downstream.stat_HETHOM_ratio,
    'INDEL_count': downstream.stat_INDEL_count,
    'sv_DUP_count': downstream.stat_sv_DUP_count,
    'sv_DEL_count': downstream.stat_sv_DEL_count,
    'sv_INS_count': downstream.stat_sv_INS_count,
    'sv_INV_count': downstream.stat_sv_INV_count,
    'sv_SWAP_count': downstream.stat_sv_SWAP_count,
    'sv_BND_count': downstream.stat_sv_BND_count,
    'trgt_genotyped_count': upstream.stat_trgt_genotyped_count,
    'trgt_uncalled_count': upstream.stat_trgt_uncalled_count
  }

  call Utilities.consolidate_stats {
    input:
      id                 = family.family_id,
      stats              = stats,
      msg_array          = flatten([flatten(upstream.msg)]),
      runtime_attributes = default_runtime_attributes
  }

  output {
    # to maintain order of samples
    Array[String] sample_ids = sample_id
    File  stats_file         = consolidate_stats.output_tsv
    File  msg_file           = consolidate_stats.messages

    # bam stats
    Array[File]   bam_statistics           = downstream.bam_statistics
    Array[File]   read_length_plot         = downstream.read_length_plot
    Array[File?]  read_quality_plot        = downstream.read_quality_plot
    Array[File]   mapq_distribution_plot   = downstream.mapq_distribution_plot
    Array[File]   mg_distribution_plot     = downstream.mg_distribution_plot
    Array[String] stat_num_reads           = downstream.stat_num_reads
    Array[String] stat_read_length_mean    = downstream.stat_read_length_mean
    Array[String] stat_read_length_median  = downstream.stat_read_length_median
    Array[String] stat_read_quality_mean   = downstream.stat_read_quality_mean
    Array[String] stat_read_quality_median = downstream.stat_read_quality_median
    Array[String] stat_mapped_read_count   = downstream.stat_mapped_read_count
    Array[String] stat_mapped_percent      = downstream.stat_mapped_percent

    # merged, haplotagged alignments
    Array[File]   merged_haplotagged_bam       = downstream.merged_haplotagged_bam
    Array[File]   merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index

    # mosdepth outputs
    Array[File]   mosdepth_summary                 = upstream.mosdepth_summary
    Array[File]   mosdepth_region_bed              = upstream.mosdepth_region_bed
    Array[File]   mosdepth_region_bed_index        = upstream.mosdepth_region_bed_index
    Array[File]   mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    Array[String] stat_mean_depth                  = upstream.stat_mean_depth
    Array[String] inferred_sex                     = upstream.inferred_sex

    # phasing stats
    Array[File]   phase_stats           = downstream.phase_stats
    Array[File]   phase_blocks          = downstream.phase_blocks
    Array[File]   phase_haplotags       = downstream.phase_haplotags
    Array[String] stat_phased_basepairs = downstream.stat_phased_basepairs
    Array[String] stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # cpg_pileup outputs
    Array[File?]  cpg_combined_bed        = downstream.cpg_combined_bed
    Array[File?]  cpg_combined_bed_index  = downstream.cpg_combined_bed_index
    Array[File?]  cpg_hap1_bed            = downstream.cpg_hap1_bed
    Array[File?]  cpg_hap1_bed_index      = downstream.cpg_hap1_bed_index
    Array[File?]  cpg_hap2_bed            = downstream.cpg_hap2_bed
    Array[File?]  cpg_hap2_bed_index      = downstream.cpg_hap2_bed_index
    Array[File?]  cpg_combined_bw         = downstream.cpg_combined_bw
    Array[File?]  cpg_hap1_bw             = downstream.cpg_hap1_bw
    Array[File?]  cpg_hap2_bw             = downstream.cpg_hap2_bw
    Array[String] stat_cpg_hap1_count     = downstream.stat_hap1_cpg_count
    Array[String] stat_cpg_hap2_count     = downstream.stat_hap2_cpg_count
    Array[String] stat_cpg_combined_count = downstream.stat_combined_cpg_count

    # sv outputs
    Array[File] phased_sv_vcf                 = downstream.phased_sv_vcf
    Array[File] phased_sv_vcf_index           = downstream.phased_sv_vcf_index
    File sv_supporting_reads                  = select_first([joint.sv_supporting_reads, upstream.sv_supporting_reads[0]])
    Array[File] sv_copynum_bedgraph           = select_first([joint.sv_copynum_bedgraph, select_all(upstream.sv_copynum_bedgraph)])
    Array[File] sv_depth_bw                   = select_first([joint.sv_depth_bw, select_all(upstream.sv_depth_bw)])
    Array[File] sv_gc_bias_corrected_depth_bw = select_first([joint.sv_gc_bias_corrected_depth_bw, select_all(upstream.sv_gc_bias_corrected_depth_bw)])
    Array[File] sv_maf_bw                     = select_first([joint.sv_maf_bw, select_all(upstream.sv_maf_bw)])

    # sv stats
    Array[String] stat_sv_DUP_count  = downstream.stat_sv_DUP_count
    Array[String] stat_sv_DEL_count  = downstream.stat_sv_DEL_count
    Array[String] stat_sv_INS_count  = downstream.stat_sv_INS_count
    Array[String] stat_sv_INV_count  = downstream.stat_sv_INV_count
    Array[String] stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    Array[String] stat_sv_BND_count  = downstream.stat_sv_BND_count

    # small variant outputs
    Array[File] phased_small_variant_vcf       = downstream.phased_small_variant_vcf
    Array[File] phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    Array[File] small_variant_gvcf             = upstream.small_variant_gvcf
    Array[File] small_variant_gvcf_index       = upstream.small_variant_gvcf_index

    # small variant stats
    Array[File]   small_variant_stats             = downstream.small_variant_stats
    Array[File]   bcftools_roh_out                = downstream.bcftools_roh_out
    Array[File]   bcftools_roh_bed                = downstream.bcftools_roh_bed
    Array[String] stat_small_variant_SNV_count    = downstream.stat_SNV_count
    Array[String] stat_small_variant_INDEL_count  = downstream.stat_INDEL_count
    Array[String] stat_small_variant_TSTV_ratio   = downstream.stat_TSTV_ratio
    Array[String] stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    Array[File]   snv_distribution_plot           = downstream.snv_distribution_plot
    Array[File]   indel_distribution_plot         = downstream.indel_distribution_plot

    # trgt outputs
    Array[File]   phased_trgt_vcf           = downstream.phased_trgt_vcf
    Array[File]   phased_trgt_vcf_index     = downstream.phased_trgt_vcf_index
    Array[File]   trgt_spanning_reads       = upstream.trgt_spanning_reads
    Array[File]   trgt_spanning_reads_index = upstream.trgt_spanning_reads_index
    Array[File]   trgt_coverage_dropouts    = downstream.trgt_coverage_dropouts
    Array[String] stat_trgt_genotyped_count = upstream.stat_trgt_genotyped_count
    Array[String] stat_trgt_uncalled_count  = upstream.stat_trgt_uncalled_count

    # paraphase outputs
    Array[File?] paraphase_output_json         = upstream.paraphase_output_json
    Array[File?] paraphase_realigned_bam       = upstream.paraphase_realigned_bam
    Array[File?] paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    Array[File?] paraphase_vcfs                = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    Array[File] mitorsaw_vcf       = upstream.mitorsaw_vcf
    Array[File] mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    Array[File] mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # PGx outputs
    Array[File]  pbstarphase_json        = downstream.pbstarphase_json
    Array[File?] pharmcat_match_json     = downstream.pharmcat_match_json
    Array[File?] pharmcat_phenotype_json = downstream.pharmcat_phenotype_json
    Array[File?] pharmcat_report_html    = downstream.pharmcat_report_html
    Array[File?] pharmcat_report_json    = downstream.pharmcat_report_json

    # joint call outputs
    File? joint_small_variants_vcf       = merge_small_variant_vcfs.merged_vcf
    File? joint_small_variants_vcf_index = merge_small_variant_vcfs.merged_vcf_index
    File? joint_sv_vcf                   = merge_sv_vcfs.merged_vcf
    File? joint_sv_vcf_index             = merge_sv_vcfs.merged_vcf_index
    File? joint_trgt_vcf                 = trgt_merge.merged_vcf
    File? joint_trgt_vcf_index           = trgt_merge.merged_vcf_index

    # tertiary analysis outputs
    File? tertiary_small_variant_filtered_vcf           = tertiary_analysis.small_variant_filtered_vcf
    File? tertiary_small_variant_filtered_vcf_index     = tertiary_analysis.small_variant_filtered_vcf_index
    File? tertiary_small_variant_filtered_tsv           = tertiary_analysis.small_variant_filtered_tsv
    File? tertiary_small_variant_compound_het_vcf       = tertiary_analysis.small_variant_compound_het_vcf
    File? tertiary_small_variant_compound_het_vcf_index = tertiary_analysis.small_variant_compound_het_vcf_index
    File? tertiary_small_variant_compound_het_tsv       = tertiary_analysis.small_variant_compound_het_tsv
    File? tertiary_sv_filtered_vcf                      = tertiary_analysis.sv_filtered_vcf
    File? tertiary_sv_filtered_vcf_index                = tertiary_analysis.sv_filtered_vcf_index
    File? tertiary_sv_filtered_tsv                      = tertiary_analysis.sv_filtered_tsv

    # qc messages
    Array[String] msg = flatten(
      [
        flatten(upstream.msg)
      ]
    )

    # workflow metadata
    String workflow_name    = "humanwgs_family"
    String workflow_version = "v3.0.1" + if defined(debug_version) then "~{"-" + debug_version}" else ""
  }
}