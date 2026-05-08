version 1.0

import "downstream/downstream.wdl" as Downstream
import "humanwgs_structs.wdl"
import "joint/joint.wdl" as Joint
import "process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "tertiary/tertiary.wdl" as TertiaryAnalysis
import "upstream/upstream.wdl" as Upstream
import "wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration

workflow humanwgs_family {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline, with joint calling for related samples."
    outputs: {
      sample_ids: {
        description: "Sample IDs"
      },
      stats_file: {
        description: "Table of summary statistics"
      },
      msg_file: {
        description: "File containing messages from the workflow"
      },
      bam_statistics: {
        description: "BAM statistics"
      },
      read_length_plot: {
        description: "Distribution of read lengths"
      },
      read_quality_plot: {
        description: "Distribution of read qualities"
      },
      mapq_distribution_plot: {
        description: "Distribution of mapping quality per alignment"
      },
      mg_distribution_plot: {
        description: "Distribution of gap-compressed identity per alignment"
      },
      stat_read_count: {
        description: "Number of reads"
      },
      stat_read_length_mean: {
        description: "Mean read length"
      },
      stat_read_length_median: {
        description: "Median read length"
      },
      stat_read_length_n50: {
        description: "Read length N50"
      },
      stat_read_quality_mean: {
        description: "Mean read quality"
      },
      stat_read_quality_median: {
        description: "Median read quality"
      },
      stat_mapped_read_count: {
        description: "Number of reads mapped to reference"
      },
      stat_mapped_read_percent: {
        description: "Percent of reads mapped to reference"
      },
      stat_gap_compressed_identity_mean: {
        description: "Mean gap-compressed identity"
      },
      stat_gap_compressed_identity_median: {
        description: "Median gap-compressed identity"
      },
      merged_haplotagged_bam: {
        description: "Merged, haplotagged alignments"
      },
      merged_haplotagged_bam_index: {
        description: "Index for merged, haplotagged alignments"
      },
      mosdepth_summary: {
        description: "Summary of aligned read depth"
      },
      mosdepth_region_bed: {
        description: "Median aligned read depth by 500bp windows"
      },
      mosdepth_region_bed_index: {
        description: "Index for median aligned read depth by 500bp windows"
      },
      mosdepth_depth_distribution_plot: {
        description: "Distribution of aligned read depth"
      },
      stat_depth_mean: {
        description: "Mean depth"
      },
      inferred_sex: {
        description: "Inferred sex"
      },
      phase_stats: {
        description: "Phasing statistics"
      },
      phase_blocks: {
        description: "Phase blocks"
      },
      phase_haplotags: {
        description: "Per-read phase assignment"
      },
      stat_phased_basepairs: {
        description: "Number of basepairs within phase blocks"
      },
      stat_phase_block_ng50: {
        description: "Phase block NG50"
      },
      cpg_combined_bed: {
        description: "5mCpG combined BED"
      },
      cpg_combined_bed_index: {
        description: "Index for 5mCpG combined BED"
      },
      cpg_hap1_bed: {
        description: "5mCpG haplotype 1 BED"
      },
      cpg_hap1_bed_index: {
        description: "Index for 5mCpG haplotype 1 BED"
      },
      cpg_hap2_bed: {
        description: "5mCpG haplotype 2 BED"
      },
      cpg_hap2_bed_index: {
        description: "Index for 5mCpG haplotype 2 BED"
      },
      cpg_combined_bw: {
        description: "5mCpG combined BigWig"
      },
      cpg_hap1_bw: {
        description: "5mCpG haplotype 1 BigWig"
      },
      cpg_hap2_bw: {
        description: "5mCpG haplotype 2 BigWig"
      },
      stat_cpg_hap1_count: {
        description: "Number of scored reference 5mCpGs in haplotype 1"
      },
      stat_cpg_hap2_count: {
        description: "Number of scored reference 5mCpGs in haplotype 2"
      },
      stat_cpg_combined_count: {
        description: "Number of scored reference 5mCpGs combined"
      },
      methbat_profile: {
        description: "MethBat 5mCpG profile"
      },
      stat_methbat_methylated_count: {
        description: "Number of profiled regions labeled as methylated"
      },
      stat_methbat_unmethylated_count: {
        description: "Number of profiled regions labeled as unmethylated"
      },
      stat_methbat_asm_count: {
        description: "Number of profiled regions labeled as having allele-specific methylation"
      },
      phased_sv_vcf: {
        description: "Phased structural variant VCF"
      },
      phased_sv_vcf_index: {
        description: "Index for phased structural variant VCF"
      },
      sv_supporting_reads: {
        description: "Supporting reads for structural variants"
      },
      sv_copynum_bedgraph: {
        description: "CNV copy number BEDGraph"
      },
      sv_depth_bw: {
        description: "CNV depth BigWig"
      },
      sv_gc_bias_corrected_depth_bw: {
        description: "CNV GC-bias corrected depth BigWig"
      },
      sv_maf_bw: {
        description: "CNV MAF BigWig"
      },
      sv_copynum_summary: {
        description: "CNV copy number summary JSON"
      },
      stat_sv_DUP_count: {
        description: "Number of DUP structural variants"
      },
      stat_sv_DEL_count: {
        description: "Number of DEL structural variants"
      },
      stat_sv_INS_count: {
        description: "Number of INS structural variants"
      },
      stat_sv_INV_count: {
        description: "Number of INV structural variants"
      },
      stat_sv_SWAP_count: {
        description: "Number of structural variant sequence swap events"
      },
      stat_sv_BND_count: {
        description: "Number of BND structural variants"
      },
      phased_small_variant_vcf: {
        description: "Phased small variant VCF"
      },
      phased_small_variant_vcf_index: {
        description: "Index for phased small variant VCF"
      },
      small_variant_gvcf: {
        description: "Small variant GVCF"
      },
      small_variant_gvcf_index: {
        description: "Index for small variant GVCF"
      },
      small_variant_stats: {
        description: "Small variant statistics"
      },
      bcftools_roh_out: {
        description: "Regions of homozygosity"
      },
      bcftools_roh_bed: {
        description: "Regions of homozygosity BED"
      },
      stat_small_variant_SNV_count: {
        description: "Number of SNVs"
      },
      stat_small_variant_INDEL_count: {
        description: "Number of INDELs"
      },
      stat_small_variant_TSTV_ratio: {
        description: "Ts/Tv ratio"
      },
      stat_small_variant_HETHOM_ratio: {
        description: "Het/Hom ratio for SNVs"
      },
      snv_distribution_plot: {
        description: "Distribution of SNVs by REF, ALT"
      },
      indel_distribution_plot: {
        description: "Distribution of indels by size"
      },
      phased_trgt_vcf: {
        description: "Phased TRGT VCF"
      },
      phased_trgt_vcf_index: {
        description: "Index for phased TRGT VCF"
      },
      trgt_spanning_reads: {
        description: "Aligned TRGT spanning reads"
      },
      trgt_spanning_reads_index: {
        description: "Index for aligned TRGT spanning reads"
      },
      trgt_coverage_dropouts: {
        description: "TRGT regions with coverage dropouts"
      },
      stat_trgt_genotyped_count: {
        description: "Number of sites genotyped by TRGT"
      },
      stat_trgt_uncalled_count: {
        description: "Number of sites ungenotyped by TRGT"
      },
      paraphase_summary: {
        description: "Paraphase summary"
      },
      paraphase_realigned_bam: {
        description: "BAM file of reads realigned by Paraphase"
      },
      paraphase_realigned_bam_index: {
        description: "Index for BAM file of reads realigned by Paraphase"
      },
      paraphase_vcfs: {
        description: "Paraphase VCFs"
      },
      mitorsaw_vcf: {
        description: "Mitochondrial variant VCF"
      },
      mitorsaw_vcf_index: {
        description: "Index for mitochondrial variant VCF"
      },
      mitorsaw_hap_stats: {
        description: "Mitochondrial haplotype statistics"
      },
      pbstarphase_summary: {
        description: "StarPhase summary"
      },
      pharmcat_match_json: {
        description: "PharmCAT match JSON"
      },
      pharmcat_phenotype_json: {
        description: "PharmCAT phenotype JSON"
      },
      pharmcat_report_html: {
        description: "PharmCAT report HTML"
      },
      pharmcat_report_json: {
        description: "PharmCAT report JSON"
      },
      joint_small_variants_vcf: {
        description: "Joint-called small variant VCF"
      },
      joint_small_variants_vcf_index: {
        description: "Index for joint-called small variant VCF"
      },
      joint_sv_vcf: {
        description: "Joint-called structural variant VCF"
      },
      joint_sv_vcf_index: {
        description: "Index for joint-called structural variant VCF"
      },
      joint_trgt_vcf: {
        description: "Joint-called TRGT VCF"
      },
      joint_trgt_vcf_index: {
        description: "Index for joint-called TRGT VCF"
      },
      tertiary_small_variant_filtered_vcf: {
        description: "Filtered, annotated small variant VCF"
      },
      tertiary_small_variant_filtered_vcf_index: {
        description: "Index for filtered, annotated small variant VCF"
      },
      tertiary_small_variant_filtered_tsv: {
        description: "Filtered, annotated small variant TSV"
      },
      tertiary_small_variant_compound_het_vcf: {
        description: "Filtered, annotated compound heterozygous small variant VCF"
      },
      tertiary_small_variant_compound_het_vcf_index: {
        description: "Index for filtered, annotated compound heterozygous small variant VCF"
      },
      tertiary_small_variant_compound_het_tsv: {
        description: "Filtered, annotated compound heterozygous small variant TSV"
      },
      tertiary_sv_filtered_vcf: {
        description: "Filtered, annotated structural variant VCF"
      },
      tertiary_sv_filtered_vcf_index: {
        description: "Index for filtered, annotated structural variant VCF"
      },
      tertiary_sv_filtered_tsv: {
        description: "Filtered, annotated structural variant TSV"
      },
      msg: {
        description: "Messages from the workflow"
      },
      workflow_name: {
        description: "Workflow name"
      },
      workflow_version: {
        description: "Workflow version"
      }
    }
  }

  parameter_meta {
    family: {
      description: "Family struct describing samples, relationships, and unaligned BAM paths"
    }
    phenotypes: {
      description: "Comma-delimited list of HPO terms for phenotypes",
      external_help: "https://hpo.jax.org"
    }
    ref_map_file: {
      description: "TSV containing reference genome file paths; must match backend"
    }
    tertiary_map_file: {
      description: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    max_reads_per_alignment_chunk: {
      description: "Maximum reads per alignment chunk"
    }
    pharmcat_min_coverage: {
      description: "Minimum coverage for PharmCAT"
    }
    glnexus_mem_gb: {
      description: "Override GLnexus memory request (GB)"
    }
    use_gpu: {
      description: "Use GPU when possible"
    }
    use_parabricks_deepvariant: {
      description: "Use Parabricks DeepVariant for small variant calling when GPU is enabled"
    }
    backend: {
      description: "Backend where the workflow will be executed",
      choices: [
        "GCP",
        "Azure",
        "AWS-HealthOmics",
        "HPC"
      ]
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      description: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      description: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      description: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      description: "Debug version for testing purposes"
    }
  }

  input {
    Family family
    String phenotypes = "HP:0000001"
    File ref_map_file
    File? tertiary_map_file
    Int max_reads_per_alignment_chunk = 500000
    Int pharmcat_min_coverage = 10
    Int glnexus_mem_gb = 60
    Boolean use_gpu = false
    Boolean use_parabricks_deepvariant = false

    # Backend configuration
    String backend
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry
    Boolean preemptible = true
    String? debug_version
  }

  call BackendConfiguration.backend_configuration { input:
    backend = backend,
    zones = zones,
    cpuPlatform = cpuPlatform,
    gpuType = gpuType,
    container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible
    then backend_configuration.spot_runtime_attributes
    else backend_configuration.on_demand_runtime_attributes

  #@ except: DeclarationName
  Map[String, String] ref_map = read_map(ref_map_file)

  call ProcessTrgtCatalog.process_trgt_catalog { input:
    trgt_catalog = ref_map["trgt_tandem_repeat_bed"],  # !FileCoercion
    ref_fasta = ref_map["fasta"],  # !FileCoercion
    ref_index = ref_map["fasta_index"],  # !FileCoercion
    default_runtime_attributes = default_runtime_attributes
  }

  Boolean single_sample = length(family.samples) == 1

  Map[String, String] pedigree_sex = {
    "MALE": "1",
    "FEMALE": "2",
    "": "."
  }

  scatter (sample in family.samples) {
    String sample_id = sample.sample_id
    #@ except: UnusedInput
    Boolean is_trio_kid = defined(sample.father_id) && defined(sample.mother_id)  # !UnusedDeclaration
    #@ except: UnusedInput
    Boolean is_duo_kid = defined(sample.father_id) != defined(sample.mother_id)  # !UnusedDeclaration

    call Upstream.upstream { input:
      sample_id = sample.sample_id,
      sex = sample.sex,
      hifi_reads = sample.hifi_reads,
      fail_reads = sample.fail_reads,
      fail_reads_bed = process_trgt_catalog.include_fail_reads_bed,
      fail_reads_bait_index = process_trgt_catalog.fail_reads_bait_index,
      ref_map_file = ref_map_file,
      max_reads_per_alignment_chunk = max_reads_per_alignment_chunk,
      single_sample = single_sample,
      use_gpu = use_gpu,
      use_parabricks_deepvariant = use_parabricks_deepvariant,
      default_runtime_attributes = default_runtime_attributes
    }

    # write sample metadata similar to pedigree format
    # family_id, sample_id, father_id, mother_id, sex, affected
    Array[String] sample_metadata = [
      family.family_id,
      sample.sample_id,
      select_first([
        sample.father_id,
        "."
      ]),
      select_first([
        sample.mother_id,
        "."
      ]),
      pedigree_sex[upstream.inferred_sex],
      if sample.affected
        then "2"
        else "1"
    ]
  }

  if (!single_sample) {
    call Joint.joint { input:
      family_id = family.family_id,
      sample_ids = sample_id,
      gvcfs = select_all(upstream.small_variant_gvcf),
      gvcf_indices = select_all(upstream.small_variant_gvcf_index),
      discover_tars = upstream.discover_tar,
      aligned_bams = upstream.aligned_hifi_reads,
      aligned_bam_indices = upstream.aligned_hifi_reads_index,
      ref_map_file = ref_map_file,
      glnexus_mem_gb = glnexus_mem_gb,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  scatter (sample_index in range(length(family.samples))) {
    call Downstream.downstream { input:
      sample_id = sample_id[sample_index],
      sex = upstream.inferred_sex[sample_index],
      aligned_hifi_reads = upstream.aligned_hifi_reads[sample_index],
      aligned_hifi_reads_index = upstream.aligned_hifi_reads_index[sample_index],
      aligned_fail_reads = upstream.aligned_fail_reads[sample_index],
      aligned_fail_reads_index = upstream.aligned_fail_reads_index[sample_index],
      trgt_catalog = process_trgt_catalog.full_catalog,
      small_variant_vcf = select_first([
        joint.split_joint_small_variant_vcfs,
        upstream.small_variant_vcf
      ])[sample_index],
      small_variant_vcf_index = select_first([
        joint.split_joint_small_variant_vcf_indices,
        upstream.small_variant_vcf_index
      ])[sample_index],
      sv_vcf = select_first([
        joint.split_joint_structural_variant_vcfs,
        select_all(upstream.sv_vcf)
      ])[sample_index],
      sv_vcf_index = select_first([
        joint.split_joint_structural_variant_vcf_indices,
        select_all(upstream.sv_vcf_index)
      ])[sample_index],
      pharmcat_min_coverage = pharmcat_min_coverage,
      ref_map_file = ref_map_file,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  if (!single_sample) {
    call Bcftools.bcftools_merge as merge_small_variant_vcfs { input:
      vcfs = downstream.phased_small_variant_vcf,
      vcf_indices = downstream.phased_small_variant_vcf_index,
      out_prefix = "~{family.family_id}.joint.~{ref_map["name"]}.small_variants.phased",
      runtime_attributes = default_runtime_attributes
    }

    call Bcftools.bcftools_merge as merge_sv_vcfs { input:
      vcfs = downstream.phased_sv_vcf,
      vcf_indices = downstream.phased_sv_vcf_index,
      out_prefix = "~{family.family_id}.joint.~{ref_map["name"]}.structural_variants.phased",
      runtime_attributes = default_runtime_attributes
    }

    call Trgt.trgt_merge { input:
      vcfs = downstream.trgt_vcf,
      vcf_indices = downstream.trgt_vcf_index,
      ref_fasta = ref_map["fasta"],  # !FileCoercion
      ref_index = ref_map["fasta_index"],  # !FileCoercion
      out_prefix = "~{family.family_id}.merged.~{ref_map["name"]}.trgt",
      runtime_attributes = default_runtime_attributes
    }
  }

  if (defined(tertiary_map_file)) {
    call TertiaryAnalysis.tertiary_analysis { input:
      sample_metadata = sample_metadata,
      phenotypes = phenotypes,
      is_trio_kid = is_trio_kid,
      is_duo_kid = is_duo_kid,
      small_variant_vcf = select_first([
        merge_small_variant_vcfs.merged_vcf,
        downstream.phased_small_variant_vcf[0]
      ]),
      small_variant_vcf_index = select_first([
        merge_small_variant_vcfs.merged_vcf_index,
        downstream.phased_small_variant_vcf_index[0]
      ]),
      sv_vcf = select_first([
        merge_sv_vcfs.merged_vcf,
        downstream.phased_sv_vcf[0]
      ]),
      sv_vcf_index = select_first([
        merge_sv_vcfs.merged_vcf_index,
        downstream.phased_sv_vcf_index[0]
      ]),
      ref_map_file = ref_map_file,
      tertiary_map_file = select_first([
        tertiary_map_file
      ]),
      default_runtime_attributes = default_runtime_attributes
    }
  }

  Array[Array[String]] stats = [
    flatten([
      [
        "sample_id"
      ],
      sample_id
    ]),
    flatten([
      [
        "read_count"
      ],
      downstream.stat_read_count
    ]),
    flatten([
      [
        "read_length_mean"
      ],
      downstream.stat_read_length_mean
    ]),
    flatten([
      [
        "read_length_median"
      ],
      downstream.stat_read_length_median
    ]),
    flatten([
      [
        "read_length_n50"
      ],
      downstream.stat_read_length_n50
    ]),
    flatten([
      [
        "read_quality_mean"
      ],
      downstream.stat_read_quality_mean
    ]),
    flatten([
      [
        "read_quality_median"
      ],
      downstream.stat_read_quality_median
    ]),
    flatten([
      [
        "mapped_read_count"
      ],
      downstream.stat_mapped_read_count
    ]),
    flatten([
      [
        "mapped_read_percent"
      ],
      downstream.stat_mapped_read_percent
    ]),
    flatten([
      [
        "gap_compressed_identity_mean"
      ],
      downstream.stat_gap_compressed_identity_mean
    ]),
    flatten([
      [
        "gap_compressed_identity_median"
      ],
      downstream.stat_gap_compressed_identity_median
    ]),
    flatten([
      [
        "depth_mean"
      ],
      upstream.stat_depth_mean
    ]),
    flatten([
      [
        "inferred_sex"
      ],
      upstream.inferred_sex
    ]),
    flatten([
      [
        "stat_phased_basepairs"
      ],
      downstream.stat_phased_basepairs
    ]),
    flatten([
      [
        "phase_block_ng50"
      ],
      downstream.stat_phase_block_ng50
    ]),
    flatten([
      [
        "cpg_combined_count"
      ],
      downstream.stat_combined_cpg_count
    ]),
    flatten([
      [
        "cpg_hap1_count"
      ],
      downstream.stat_hap1_cpg_count
    ]),
    flatten([
      [
        "cpg_hap2_count"
      ],
      downstream.stat_hap2_cpg_count
    ]),
    flatten([
      [
        "methbat_methylated_count"
      ],
      downstream.stat_methbat_methylated_count
    ]),
    flatten([
      [
        "methbat_unmethylated_count"
      ],
      downstream.stat_methbat_unmethylated_count
    ]),
    flatten([
      [
        "methbat_asm_count"
      ],
      downstream.stat_methbat_asm_count
    ]),
    flatten([
      [
        "SNV_count"
      ],
      downstream.stat_SNV_count
    ]),
    flatten([
      [
        "TSTV_ratio"
      ],
      downstream.stat_TSTV_ratio
    ]),
    flatten([
      [
        "HETHOM_ratio"
      ],
      downstream.stat_HETHOM_ratio
    ]),
    flatten([
      [
        "INDEL_count"
      ],
      downstream.stat_INDEL_count
    ]),
    flatten([
      [
        "sv_DUP_count"
      ],
      downstream.stat_sv_DUP_count
    ]),
    flatten([
      [
        "sv_DEL_count"
      ],
      downstream.stat_sv_DEL_count
    ]),
    flatten([
      [
        "sv_INS_count"
      ],
      downstream.stat_sv_INS_count
    ]),
    flatten([
      [
        "sv_INV_count"
      ],
      downstream.stat_sv_INV_count
    ]),
    flatten([
      [
        "sv_SWAP_count"
      ],
      downstream.stat_sv_SWAP_count
    ]),
    flatten([
      [
        "sv_BND_count"
      ],
      downstream.stat_sv_BND_count
    ]),
    flatten([
      [
        "trgt_genotyped_count"
      ],
      downstream.stat_trgt_genotyped_count
    ]),
    flatten([
      [
        "trgt_uncalled_count"
      ],
      downstream.stat_trgt_uncalled_count
    ])
  ]

  call Utilities.consolidate_stats { input:
    out_prefix = family.family_id,
    stats = stats,
    msg_array = flatten([
      process_trgt_catalog.msg,
      flatten(upstream.msg)
    ]),
    runtime_attributes = default_runtime_attributes
  }

  output {
    # to maintain order of samples
    Array[String] sample_ids = sample_id
    File stats_file = consolidate_stats.stats_tsv
    File msg_file = consolidate_stats.messages

    # bam stats
    Array[File] bam_statistics = downstream.bam_statistics
    Array[File] read_length_plot = downstream.read_length_plot
    Array[File?] read_quality_plot = downstream.read_quality_plot
    Array[File] mapq_distribution_plot = downstream.mapq_distribution_plot
    Array[File] mg_distribution_plot = downstream.mg_distribution_plot
    Array[String] stat_read_count = downstream.stat_read_count
    Array[String] stat_read_length_mean = downstream.stat_read_length_mean
    Array[String] stat_read_length_median = downstream.stat_read_length_median
    Array[String] stat_read_length_n50 = downstream.stat_read_length_n50
    Array[String] stat_read_quality_mean = downstream.stat_read_quality_mean
    Array[String] stat_read_quality_median = downstream.stat_read_quality_median
    Array[String] stat_mapped_read_count = downstream.stat_mapped_read_count
    Array[String] stat_mapped_read_percent = downstream.stat_mapped_read_percent
    Array[String] stat_gap_compressed_identity_mean = downstream.stat_gap_compressed_identity_mean
    Array[String] stat_gap_compressed_identity_median = downstream.stat_gap_compressed_identity_median

    # merged, haplotagged alignments
    Array[File] merged_haplotagged_bam = downstream.merged_haplotagged_bam
    Array[File] merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index

    # mosdepth outputs
    Array[File] mosdepth_summary = upstream.mosdepth_summary
    Array[File] mosdepth_region_bed = upstream.mosdepth_region_bed
    Array[File] mosdepth_region_bed_index = upstream.mosdepth_region_bed_index
    Array[File] mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    Array[String] stat_depth_mean = upstream.stat_depth_mean
    Array[String] inferred_sex = upstream.inferred_sex

    # phasing stats
    Array[File] phase_stats = downstream.phase_stats
    Array[File] phase_blocks = downstream.phase_blocks
    Array[File] phase_haplotags = downstream.phase_haplotags
    Array[String] stat_phased_basepairs = downstream.stat_phased_basepairs
    Array[String] stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # methylation outputs and profile
    Array[File?] cpg_combined_bed = downstream.cpg_combined_bed
    Array[File?] cpg_combined_bed_index = downstream.cpg_combined_bed_index
    Array[File?] cpg_hap1_bed = downstream.cpg_hap1_bed
    Array[File?] cpg_hap1_bed_index = downstream.cpg_hap1_bed_index
    Array[File?] cpg_hap2_bed = downstream.cpg_hap2_bed
    Array[File?] cpg_hap2_bed_index = downstream.cpg_hap2_bed_index
    Array[File?] cpg_combined_bw = downstream.cpg_combined_bw
    Array[File?] cpg_hap1_bw = downstream.cpg_hap1_bw
    Array[File?] cpg_hap2_bw = downstream.cpg_hap2_bw
    Array[String] stat_cpg_hap1_count = downstream.stat_hap1_cpg_count
    Array[String] stat_cpg_hap2_count = downstream.stat_hap2_cpg_count
    Array[String] stat_cpg_combined_count = downstream.stat_combined_cpg_count
    Array[File?] methbat_profile = downstream.methbat_profile
    Array[String] stat_methbat_methylated_count = downstream.stat_methbat_methylated_count
    Array[String] stat_methbat_unmethylated_count = downstream.stat_methbat_unmethylated_count
    Array[String] stat_methbat_asm_count = downstream.stat_methbat_asm_count

    # sv outputs
    Array[File] phased_sv_vcf = downstream.phased_sv_vcf
    Array[File] phased_sv_vcf_index = downstream.phased_sv_vcf_index
    File sv_supporting_reads = select_first([
      joint.sv_supporting_reads,
      upstream.sv_supporting_reads[0]
    ])
    Array[File] sv_copynum_bedgraph = select_first([
      joint.sv_copynum_bedgraph,
      select_all(upstream.sv_copynum_bedgraph)
    ])
    Array[File] sv_depth_bw = select_first([
      joint.sv_depth_bw,
      select_all(upstream.sv_depth_bw)
    ])
    Array[File] sv_gc_bias_corrected_depth_bw = select_first([
      joint.sv_gc_bias_corrected_depth_bw,
      select_all(upstream.sv_gc_bias_corrected_depth_bw)
    ])
    Array[File] sv_maf_bw = select_first([
      joint.sv_maf_bw,
      select_all(upstream.sv_maf_bw)
    ])
    Array[File] sv_copynum_summary = select_first([
      joint.sv_copynum_summary,
      select_all(upstream.sv_copynum_summary)
    ])

    # sv stats
    Array[String] stat_sv_DUP_count = downstream.stat_sv_DUP_count
    Array[String] stat_sv_DEL_count = downstream.stat_sv_DEL_count
    Array[String] stat_sv_INS_count = downstream.stat_sv_INS_count
    Array[String] stat_sv_INV_count = downstream.stat_sv_INV_count
    Array[String] stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    Array[String] stat_sv_BND_count = downstream.stat_sv_BND_count

    # small variant outputs
    Array[File] phased_small_variant_vcf = downstream.phased_small_variant_vcf
    Array[File] phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    Array[File?] small_variant_gvcf = upstream.small_variant_gvcf
    Array[File?] small_variant_gvcf_index = upstream.small_variant_gvcf_index

    # small variant stats
    Array[File] small_variant_stats = downstream.small_variant_stats
    Array[File] bcftools_roh_out = downstream.bcftools_roh_out
    Array[File] bcftools_roh_bed = downstream.bcftools_roh_bed
    Array[String] stat_small_variant_SNV_count = downstream.stat_SNV_count
    Array[String] stat_small_variant_INDEL_count = downstream.stat_INDEL_count
    Array[String] stat_small_variant_TSTV_ratio = downstream.stat_TSTV_ratio
    Array[String] stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    Array[File] snv_distribution_plot = downstream.snv_distribution_plot
    Array[File] indel_distribution_plot = downstream.indel_distribution_plot

    # trgt outputs
    Array[File] phased_trgt_vcf = downstream.trgt_vcf
    Array[File] phased_trgt_vcf_index = downstream.trgt_vcf_index
    Array[File] trgt_spanning_reads = downstream.trgt_spanning_reads
    Array[File] trgt_spanning_reads_index = downstream.trgt_spanning_reads_index
    Array[File] trgt_coverage_dropouts = downstream.trgt_coverage_dropouts
    Array[String] stat_trgt_genotyped_count = downstream.stat_trgt_genotyped_count
    Array[String] stat_trgt_uncalled_count = downstream.stat_trgt_uncalled_count

    # paraphase outputs
    Array[File?] paraphase_summary = upstream.paraphase_output_json
    Array[File?] paraphase_realigned_bam = upstream.paraphase_realigned_bam
    Array[File?] paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    Array[File?] paraphase_vcfs = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    Array[File] mitorsaw_vcf = upstream.mitorsaw_vcf
    Array[File] mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    Array[File] mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # PGx outputs
    Array[File] pbstarphase_summary = downstream.pbstarphase_json
    Array[File?] pharmcat_match_json = downstream.pharmcat_match_json
    Array[File?] pharmcat_phenotype_json = downstream.pharmcat_phenotype_json
    Array[File?] pharmcat_report_html = downstream.pharmcat_report_html
    Array[File?] pharmcat_report_json = downstream.pharmcat_report_json

    # joint call outputs
    File? joint_small_variants_vcf = merge_small_variant_vcfs.merged_vcf
    File? joint_small_variants_vcf_index = merge_small_variant_vcfs.merged_vcf_index
    File? joint_sv_vcf = merge_sv_vcfs.merged_vcf
    File? joint_sv_vcf_index = merge_sv_vcfs.merged_vcf_index
    File? joint_trgt_vcf = trgt_merge.merged_vcf
    File? joint_trgt_vcf_index = trgt_merge.merged_vcf_index

    # tertiary analysis outputs
    File? tertiary_small_variant_filtered_vcf = tertiary_analysis.small_variant_filtered_vcf
    File? tertiary_small_variant_filtered_vcf_index = tertiary_analysis.small_variant_filtered_vcf_index
    File? tertiary_small_variant_filtered_tsv = tertiary_analysis.small_variant_filtered_tsv
    File? tertiary_small_variant_compound_het_vcf = tertiary_analysis.small_variant_compound_het_vcf
    File? tertiary_small_variant_compound_het_vcf_index = tertiary_analysis.small_variant_compound_het_vcf_index
    File? tertiary_small_variant_compound_het_tsv = tertiary_analysis.small_variant_compound_het_tsv
    File? tertiary_sv_filtered_vcf = tertiary_analysis.sv_filtered_vcf
    File? tertiary_sv_filtered_vcf_index = tertiary_analysis.sv_filtered_vcf_index
    File? tertiary_sv_filtered_tsv = tertiary_analysis.sv_filtered_tsv

    # qc messages
    Array[String] msg = flatten([
      process_trgt_catalog.msg,
      flatten(upstream.msg),
      flatten(downstream.msg)
    ])

    # workflow metadata
    String workflow_name = "humanwgs_family"
    String workflow_version = "v3.3.1" + if defined(debug_version)
      then "~{"-" + debug_version}"
      else ""
  }
}
