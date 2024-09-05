version 1.0

import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "upstream/upstream.wdl" as Upstream
import "joint/joint.wdl" as Joint
import "downstream/downstream.wdl" as Downstream
import "wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "wdl-common/wdl/tasks/write_ped_phrank.wdl" as Write_ped_phrank
import "tertiary/tertiary.wdl" as TertiaryAnalysis


workflow humanwgs_family {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline, with joint calling for related samples."
  }

  parameter_meta {
    family: {
      name: "Family struct describing samples, relationships, and unaligned BAM paths"
    }
    ref_map_file: {
      name: "TSV containing reference genome file paths; must match backend"
    }
    deepvariant_version: {
      name: "DeepVariant version"
    }
    custom_deepvariant_model_tar: {
      name: "Custom DeepVariant model tarball"
    }
    pharmcat_version: {
      name: "PharmCAT version"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    tertiary_map_file: {
      name: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    glnexus_mem_gb: {
      name: "Override GLnexus memory request (GB)"
    }
    pbsv_call_mem_gb: {
      name: "Override PBSV call memory request (GB)"
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
    gpuType: {
      name: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      name: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'"
    }
    container_namespace: {
      name: "AWS ECRs have the format REGISTRY/NAMESPACE/CONTAINER. Must be set if backend is set to 'AWS-HealthOmics'"
    }
    preemptible: {
      name: "Where possible, run tasks preemptibly"
    }
  }

  input {
    Family family

    File ref_map_file

    # These options are only intended for testing purposes.
    # There is no guarantee that the pipeline will work with
    # other version of DeepVariant or with custom models.
    String deepvariant_version = "1.6.1"
    File? custom_deepvariant_model_tar

    String pharmcat_version = "2.15.0"
    Int pharmcat_min_coverage = 10

    String phenotypes = "HP:0000001"
    File? tertiary_map_file

    Int? glnexus_mem_gb
    Int? pbsv_call_mem_gb

    Boolean gpu = false

    # Backend configuration
    String backend
    String? zones
    String? gpuType
    String? container_registry
    String? container_namespace

    Boolean preemptible = true
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = backend,
      zones              = zones,
      gpuType            = gpuType,
      container_registry = if defined(container_namespace) then select_first([container_registry]) + "/" + select_first([container_namespace]) else container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  Map [String, String] ref_map = read_map(ref_map_file)

  scatter (sample in family.samples) {
    String sample_id = sample.sample_id
    call Upstream.upstream {
      input:
        sample_id                    = sample.sample_id,
        sex                          = sample.sex,
        hifi_reads                   = sample.hifi_reads,
        ref_map_file                 = ref_map_file,
        deepvariant_version          = deepvariant_version,
        custom_deepvariant_model_tar = custom_deepvariant_model_tar,
        gpu                          = gpu,
        default_runtime_attributes   = default_runtime_attributes
    }
  }

  call Joint.joint {
    input:
      family_id                  = family.family_id,
      sample_ids                 = sample_id,
      gvcfs                      = upstream.small_variant_vcf,
      gvcf_indices               = upstream.small_variant_vcf_index,
      svsigs                     = flatten(upstream.svsigs),
      ref_map_file               = ref_map_file,
      glnexus_mem_gb             = glnexus_mem_gb,
      pbsv_call_mem_gb           = pbsv_call_mem_gb,
      default_runtime_attributes = default_runtime_attributes
  }

  scatter (sample_index in range(length(family.samples))) {
    call Downstream.downstream {
      input:
        sample_id                  = sample_id[sample_index],
        small_variant_vcf          = joint.split_joint_small_variant_vcfs[sample_index],
        small_variant_vcf_index    = joint.split_joint_small_variant_vcf_indices[sample_index],
        sv_vcf                     = joint.split_joint_structural_variant_vcfs[sample_index],
        sv_vcf_index               = joint.split_joint_structural_variant_vcf_indices[sample_index],
        trgt_vcf                   = upstream.trgt_vcf[sample_index],
        trgt_vcf_index             = upstream.trgt_vcf_index[sample_index],
        aligned_bam                = upstream.out_bam[sample_index],
        aligned_bam_index          = upstream.out_bam_index[sample_index],
        pharmcat_version           = pharmcat_version,
        pharmcat_min_coverage      = pharmcat_min_coverage,
        ref_map_file               = ref_map_file,
        default_runtime_attributes = default_runtime_attributes
    }
  }

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
      out_prefix         = "~{family.family_id}.~{ref_map['name']}.trgt",
      runtime_attributes = default_runtime_attributes
  }

  if (defined(tertiary_map_file)) {
    scatter (sample_index in range(length(family.samples))) {
      Map[String, String] sample_metadata = {
        "sample_id": family.samples[sample_index].sample_id,
        "sex": select_first([family.samples[sample_index].sex, upstream.inferred_sex[sample_index]]),
        "affected": family.samples[sample_index].affected,
        "father_id": select_first([family.samples[sample_index].father_id, "."]),
        "mother_id": select_first([family.samples[sample_index].mother_id, "."])
      }
    }

    call Write_ped_phrank.write_ped_phrank {
      input:
        id                 = family.family_id,
        family_json        = sample_metadata,
        phenotypes         = phenotypes,
        runtime_attributes = default_runtime_attributes
    }

    call TertiaryAnalysis.tertiary_analysis {
      input:
        pedigree                   = write_ped_phrank.pedigree,
        phrank_lookup              = write_ped_phrank.phrank_lookup,
        small_variant_vcf          = merge_small_variant_vcfs.merged_vcf,
        small_variant_vcf_index    = merge_small_variant_vcfs.merged_vcf_index,
        sv_vcf                     = merge_sv_vcfs.merged_vcf,
        sv_vcf_index               = merge_sv_vcfs.merged_vcf_index,
        ref_map_file               = ref_map_file,
        tertiary_map_file          = select_first([tertiary_map_file]),
        default_runtime_attributes = default_runtime_attributes
    }
  }

  output {
    # to maintain order of samples
    Array[String] sample_ids = sample_id

    # bam stats
    Array[File] bam_stats                  = upstream.read_length_and_quality
    Array[File] read_length_histogram      = upstream.read_length_histogram
    Array[File] read_quality_histogram     = upstream.read_quality_histogram
    Array[File] read_length_plot           = upstream.read_length_plot
    Array[File] read_quality_plot          = upstream.read_quality_plot
    Array[String] stat_num_reads           = upstream.stat_num_reads
    Array[String] stat_read_length_mean    = upstream.stat_read_length_mean
    Array[String] stat_read_length_median  = upstream.stat_read_length_median
    Array[String] stat_read_quality_mean   = upstream.stat_read_quality_mean
    Array[String] stat_read_quality_median = upstream.stat_read_quality_median

    # merged, haplotagged alignments
    Array[File] merged_haplotagged_bam       = downstream.merged_haplotagged_bam
    Array[File] merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index
    Array[String] stat_mapped_fraction       = downstream.stat_mapped_fraction

    # mosdepth outputs
    Array[File] mosdepth_summary          = upstream.mosdepth_summary
    Array[File] mosdepth_region_bed       = upstream.mosdepth_region_bed
    Array[File] mosdepth_region_bed_index = upstream.mosdepth_region_bed
    Array[String] inferred_sex            = upstream.inferred_sex
    Array[String] stat_mean_depth         = upstream.stat_mean_depth

    # phasing stats
    Array[File] phase_stats             = downstream.phase_stats
    Array[File] phase_blocks            = downstream.phase_blocks
    Array[File] phase_haplotags         = downstream.phase_haplotags
    Array[String] stat_phased_basepairs = downstream.stat_phased_basepairs
    Array[String] stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # cpg_pileup outputs
    Array[File] cpg_combined_bed          = downstream.cpg_combined_bed
    Array[File] cpg_hap1_bed              = downstream.cpg_hap1_bed
    Array[File] cpg_hap2_bed              = downstream.cpg_hap2_bed
    Array[File] cpg_combined_bw           = downstream.cpg_combined_bw
    Array[File] cpg_hap1_bw               = downstream.cpg_hap1_bw
    Array[File] cpg_hap2_bw               = downstream.cpg_hap2_bw
    Array[String] stat_hap1_cpg_count     = downstream.stat_hap1_cpg_count
    Array[String] stat_hap2_cpg_count     = downstream.stat_hap2_cpg_count
    Array[String] stat_combined_cpg_count = downstream.stat_combined_cpg_count

    # sv outputs
    Array[File] phased_sv_vcf       = downstream.phased_sv_vcf
    Array[File] phased_sv_vcf_index = downstream.phased_sv_vcf_index

    # sv stats
    Array[String] stat_sv_DUP_count = downstream.stat_sv_DUP_count
    Array[String] stat_sv_DEL_count = downstream.stat_sv_DEL_count
    Array[String] stat_sv_INS_count = downstream.stat_sv_INS_count
    Array[String] stat_sv_INV_count = downstream.stat_sv_INV_count
    Array[String] stat_sv_BND_count = downstream.stat_sv_BND_count

    # small variant outputs
    Array[File] phased_small_variant_vcf       = downstream.phased_small_variant_vcf
    Array[File] phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    Array[File] small_variant_gvcf             = upstream.small_variant_gvcf
    Array[File] small_variant_gvcf_index       = upstream.small_variant_gvcf_index

    # small variant stats
    Array[File] small_variant_stats = downstream.small_variant_stats
    Array[File] bcftools_roh_out    = downstream.bcftools_roh_out
    Array[File] bcftools_roh_bed    = downstream.bcftools_roh_bed
    Array[String] stat_SNV_count    = downstream.stat_SNV_count
    Array[String] stat_INDEL_count  = downstream.stat_INDEL_count
    Array[String] stat_TSTV_ratio   = downstream.stat_TSTV_ratio
    Array[String] stat_HETHOM_ratio = downstream.stat_HETHOM_ratio

    # trgt outputs
    Array[File] phased_trgt_vcf             = downstream.phased_trgt_vcf
    Array[File] phased_trgt_vcf_index       = downstream.phased_trgt_vcf_index
    Array[File] trgt_spanning_reads         = upstream.trgt_spanning_reads
    Array[File] trgt_spanning_reads_index   = upstream.trgt_spanning_reads_index
    Array[File] trgt_coverage_dropouts      = upstream.trgt_coverage_dropouts
    Array[String] stat_trgt_genotyped_count = upstream.stat_trgt_genotyped_count
    Array[String] stat_trgt_uncalled_count  = upstream.stat_trgt_uncalled_count

    # paraphase outputs
    Array[File] paraphase_output_json         = upstream.paraphase_output_json
    Array[File] paraphase_realigned_bam       = upstream.paraphase_realigned_bam
    Array[File] paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    Array[File?] paraphase_vcfs               = upstream.paraphase_vcfs

    # per sample ficnv outputs
    Array[File] cnv_vcf              = upstream.cnv_vcf
    Array[File] cnv_vcf_index        = upstream.cnv_vcf_index
    Array[File] cnv_copynum_bedgraph = upstream.cnv_copynum_bedgraph
    Array[File] cnv_depth_bw         = upstream.cnv_depth_bw
    Array[File] cnv_maf_bw           = upstream.cnv_maf_bw
    Array[String] stat_cnv_DUP_count = upstream.stat_cnv_DUP_count
    Array[String] stat_cnv_DEL_count = upstream.stat_cnv_DEL_count
    Array[String] stat_cnv_DUP_sum   = upstream.stat_cnv_DUP_sum
    Array[String] stat_cnv_DEL_sum   = upstream.stat_cnv_DEL_sum

    # PGx outputs
    Array[File] pbstarphase_json                   = downstream.pbstarphase_json
    Array[File?] pharmcat_missing_pgx_vcf          = downstream.pharmcat_missing_pgx_vcf
    Array[File] pharmcat_preprocessed_filtered_vcf = downstream.pharmcat_preprocessed_filtered_vcf
    Array[File] pharmcat_match_json                = downstream.pharmcat_match_json
    Array[File] pharmcat_phenotype_json            = downstream.pharmcat_phenotype_json
    Array[File] pharmcat_report_html               = downstream.pharmcat_report_html
    Array[File] pharmcat_report_json               = downstream.pharmcat_report_json

    # joint call outputs
    File joint_small_variants_vcf            = merge_small_variant_vcfs.merged_vcf
    File joint_small_variants_vcf_index      = merge_small_variant_vcfs.merged_vcf_index
    File joint_structural_variants_vcf       = merge_sv_vcfs.merged_vcf
    File joint_structural_variants_vcf_index = merge_sv_vcfs.merged_vcf_index
    File joint_trgt_vcf                      = trgt_merge.merged_vcf
    File joint_trgt_vcf_index                = trgt_merge.merged_vcf_index

    # tertiary analysis outputs
    File? filtered_small_variant_vcf           = tertiary_analysis.filtered_small_variant_vcf
    File? filtered_small_variant_vcf_index     = tertiary_analysis.filtered_small_variant_vcf_index
    File? filtered_small_variant_tsv           = tertiary_analysis.filtered_small_variant_tsv
    File? compound_het_small_variant_vcf       = tertiary_analysis.compound_het_small_variant_vcf
    File? compound_het_small_variant_vcf_index = tertiary_analysis.compound_het_small_variant_vcf_index
    File? compound_het_small_variant_tsv       = tertiary_analysis.compound_het_small_variant_tsv
    File? filtered_svpack_vcf                  = tertiary_analysis.filtered_svpack_vcf
    File? filtered_svpack_vcf_index            = tertiary_analysis.filtered_svpack_vcf_index
    File? filtered_svpack_tsv                  = tertiary_analysis.filtered_svpack_tsv
  }
}