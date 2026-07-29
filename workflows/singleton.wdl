version 1.0

import "downstream/downstream.wdl" as Downstream
import "process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "unpack_container_manifest/unpack_container_manifest.wdl" as UnpackContainerManifest
import "upstream/upstream.wdl" as Upstream
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration

workflow humanwgs_singleton {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline for individual samples."
    outputs: {
      stats_file: {
        description: "Table of summary statistics"
      },
      msg_file: {
        description: "File containing messages from the workflow"
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
      cpg_pileup_bed: {
        description: "5mCpG pileup BED"
      },
      cpg_pileup_bed_index: {
        description: "Index for 5mCpG pileup BED"
      },
      hmcpg_pileup_bed: {
        description: "5hmC pileup BED"
      },
      hmcpg_pileup_bed_index: {
        description: "Index for 5hmC pileup BED"
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
      sv_stats_plot: {
        description: "Structural variant size distribution plot"
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
      kivvi_kiv2_vcf: {
        description: "KIV2 repeat variant VCF"
      },
      kivvi_kiv2_vcf_index: {
        description: "Index for KIV2 repeat variant VCF"
      },
      kivvi_kiv2_json: {
        description: "KIV2 repeat genotype JSON"
      },
      kivvi_kiv2_realigned_bam: {
        description: "KIV2-realigned BAM"
      },
      kivvi_kiv2_realigned_bam_index: {
        description: "Index for KIV2-realigned BAM"
      },
      kivvi_kiv2_allele_plot: {
        description: "KIV2 assembled allele plot"
      },
      kivvi_d4z4_vcf: {
        description: "D4Z4 repeat variant VCF"
      },
      kivvi_d4z4_vcf_index: {
        description: "Index for D4Z4 repeat variant VCF"
      },
      kivvi_d4z4_json: {
        description: "D4Z4 repeat genotype JSON"
      },
      kivvi_d4z4_realigned_bam: {
        description: "D4Z4-realigned BAM"
      },
      kivvi_d4z4_realigned_bam_index: {
        description: "Index for D4Z4-realigned BAM"
      },
      kivvi_d4z4_allele_plot: {
        description: "D4Z4 assembled allele plot"
      },
      pbstarphase_summary: {
        description: "StarPhase summary"
      },
      pbstarphase_tsv: {
        description: "StarPhase summary in TSV format for PharmCAT"
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
    sample_id: {
      description: "Unique identifier for the sample"
    }
    hifi_reads: {
      description: "Array of paths to hifi_reads in unaligned BAM format"
    }
    fail_reads: {
      description: "Array of paths to fail_reads in unaligned BAM format"
    }
    ref_name: {
      description: "Reference genome to use for this workflow run",
      choices: [
        "GRCh38",
        "GRCh38_GIABv3"
      ]
    }
    trgt_tandem_repeat_bed_override: {
      description: "Optional BED file to override the default TRGT tandem repeat catalog"
    }
    methbat_region_tsv_override: {
      description: "Optional TSV file to override the default MethBat methylation profiling regions"
    }
    use_alignment_chunking: {
      description: "Whether to chunk BAM files for alignment. If false, all reads will be aligned in a single chunk."
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
      description: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'",
      default: "quay.io/pacbio"
    }
    preemptible: {
      description: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      description: "Debug version for testing purposes"
    }
  }

  input {
    String sample_id
    Array[File] hifi_reads
    Array[File]? fail_reads
    String ref_name = "GRCh38"
    File? trgt_tandem_repeat_bed_override
    File? methbat_region_tsv_override
    Boolean use_alignment_chunking = true
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

  # Pinned reference data container images, keyed by ref_name; add an entry
  # here when a new reference is supported (see upstream.wdl's
  # max_norm_female_chrY_depth for the sibling per-reference map).
  Map[String, String] reference_container = {
    "GRCh38": "workflow-data-container-hifi-human-wgs-wdl-grch38@sha256:5e3f23e44c1c09762838e81677f7d136367be2a4063efc608129841cef745b42",  # v4.0.0
    "GRCh38_GIABv3": "workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3@sha256:f3799c1eff1b816a95ea06ce567d3324d06e945428567077d96eca07c76ad0aa"  # v4.0.0
  }

  call UnpackContainerManifest.unpack_container_manifest { input:
    unpack_image = "~{default_runtime_attributes.container_registry}/~{reference_container[ref_name]}",
    runtime_attributes = default_runtime_attributes
  }

  File ref_fasta = unpack_container_manifest.ref_fasta
  File ref_index = unpack_container_manifest.ref_index
  Float max_norm_female_chrY_depth = unpack_container_manifest.max_norm_female_chrY_depth
  File trgt_tandem_repeat_bed = select_first([
    trgt_tandem_repeat_bed_override,
    unpack_container_manifest.trgt_tandem_repeat_bed
  ])
  File sawfish_exclude_bed = unpack_container_manifest.sawfish_exclude_bed
  File sawfish_exclude_bed_index = unpack_container_manifest.sawfish_exclude_bed_index
  File sawfish_expected_bed_male = unpack_container_manifest.sawfish_expected_bed_male
  File sawfish_expected_bed_female = unpack_container_manifest.sawfish_expected_bed_female
  File methbat_region_tsv = select_first([
    methbat_region_tsv_override,
    unpack_container_manifest.methbat_region_tsv
  ])
  String paraphase_genome_build = unpack_container_manifest.paraphase_genome_build
  Boolean run_starphase = unpack_container_manifest.run_starphase

  call ProcessTrgtCatalog.process_trgt_catalog { input:
    trgt_catalog = trgt_tandem_repeat_bed,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    default_runtime_attributes = default_runtime_attributes
  }

  call Upstream.upstream { input:
    sample_id = sample_id,
    hifi_reads = hifi_reads,
    fail_reads = fail_reads,
    fail_reads_bed = process_trgt_catalog.include_fail_reads_bed,
    fail_reads_bait_index = process_trgt_catalog.fail_reads_bait_index,
    ref_name = ref_name,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    max_norm_female_chrY_depth = max_norm_female_chrY_depth,
    paraphase_genome_build = paraphase_genome_build,
    sawfish_exclude_bed = sawfish_exclude_bed,
    sawfish_exclude_bed_index = sawfish_exclude_bed_index,
    sawfish_expected_bed_male = sawfish_expected_bed_male,
    sawfish_expected_bed_female = sawfish_expected_bed_female,
    use_alignment_chunking = use_alignment_chunking,
    single_sample = true,
    use_gpu = use_gpu,
    use_parabricks_deepvariant = use_parabricks_deepvariant,
    default_runtime_attributes = default_runtime_attributes
  }

  call Downstream.downstream { input:
    sample_id = sample_id,
    sex = upstream.inferred_sex,
    aligned_hifi_reads = upstream.aligned_hifi_reads,
    aligned_hifi_reads_index = upstream.aligned_hifi_reads_index,
    aligned_fail_reads = upstream.aligned_fail_reads,
    aligned_fail_reads_index = upstream.aligned_fail_reads_index,
    trgt_catalog = process_trgt_catalog.full_catalog,
    small_variant_vcf = upstream.small_variant_vcf,
    small_variant_vcf_index = upstream.small_variant_vcf_index,
    sv_vcf = select_first([
      upstream.sv_vcf
    ]),
    sv_vcf_index = select_first([
      upstream.sv_vcf_index
    ]),
    ref_name = ref_name,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    sawfish_expected_bed_male = sawfish_expected_bed_male,
    sawfish_expected_bed_female = sawfish_expected_bed_female,
    methbat_region_tsv = methbat_region_tsv,
    run_starphase = run_starphase,
    default_runtime_attributes = default_runtime_attributes
  }

  Array[Array[String]] stats = [
    [
      "sample_id",
      sample_id
    ],
    [
      "read_count",
      downstream.stat_read_count
    ],
    [
      "read_length_mean",
      downstream.stat_read_length_mean
    ],
    [
      "read_length_median",
      downstream.stat_read_length_median
    ],
    [
      "read_length_n50",
      downstream.stat_read_length_n50
    ],
    [
      "read_quality_mean",
      downstream.stat_read_quality_mean
    ],
    [
      "read_quality_median",
      downstream.stat_read_quality_median
    ],
    [
      "mapped_read_count",
      downstream.stat_mapped_read_count
    ],
    [
      "mapped_read_percent",
      downstream.stat_mapped_read_percent
    ],
    [
      "gap_compressed_identity_mean",
      downstream.stat_gap_compressed_identity_mean
    ],
    [
      "gap_compressed_identity_median",
      downstream.stat_gap_compressed_identity_median
    ],
    [
      "depth_mean",
      upstream.stat_depth_mean
    ],
    [
      "inferred_sex",
      upstream.inferred_sex
    ],
    [
      "stat_phased_basepairs",
      downstream.stat_phased_basepairs
    ],
    [
      "phase_block_ng50",
      downstream.stat_phase_block_ng50
    ],
    [
      "cpg_combined_count",
      downstream.stat_cpg_combined_count
    ],
    [
      "cpg_hap1_count",
      downstream.stat_cpg_hap1_count
    ],
    [
      "cpg_hap2_count",
      downstream.stat_cpg_hap2_count
    ],
    [
      "methbat_methylated_count",
      downstream.stat_methbat_methylated_count
    ],
    [
      "methbat_unmethylated_count",
      downstream.stat_methbat_unmethylated_count
    ],
    [
      "methbat_asm_count",
      downstream.stat_methbat_asm_count
    ],
    [
      "SNV_count",
      downstream.stat_SNV_count
    ],
    [
      "TSTV_ratio",
      downstream.stat_TSTV_ratio
    ],
    [
      "HETHOM_ratio",
      downstream.stat_HETHOM_ratio
    ],
    [
      "INDEL_count",
      downstream.stat_INDEL_count
    ],
    [
      "sv_DUP_count",
      downstream.stat_sv_DUP_count
    ],
    [
      "sv_DEL_count",
      downstream.stat_sv_DEL_count
    ],
    [
      "sv_INS_count",
      downstream.stat_sv_INS_count
    ],
    [
      "sv_INV_count",
      downstream.stat_sv_INV_count
    ],
    [
      "sv_SWAP_count",
      downstream.stat_sv_SWAP_count
    ],
    [
      "sv_BND_count",
      downstream.stat_sv_BND_count
    ],
    [
      "trgt_genotyped_count",
      downstream.stat_trgt_genotyped_count
    ],
    [
      "trgt_uncalled_count",
      downstream.stat_trgt_uncalled_count
    ]
  ]

  # msg outputs from tasks may contain empty strings when an empty
  # messages.txt is read back as [""] rather than []; drop them here
  scatter (m in flatten([
    process_trgt_catalog.msg,
    upstream.msg
  ])) {
    if (m != "") {
      String non_empty_stats_msg = m
    }
  }

  call Utilities.consolidate_stats { input:
    out_prefix = sample_id,
    stats = stats,
    msg_array = select_all(non_empty_stats_msg),
    runtime_attributes = default_runtime_attributes
  }

  scatter (m in flatten([
    process_trgt_catalog.msg,
    upstream.msg,
    downstream.msg
  ])) {
    if (m != "") {
      String non_empty_msg = m
    }
  }

  output {
    # consolidated stats
    File stats_file = consolidate_stats.stats_tsv
    File msg_file = consolidate_stats.messages

    # bam stats
    File read_length_plot = downstream.read_length_plot
    File read_quality_plot = downstream.read_quality_plot
    File mapq_distribution_plot = downstream.mapq_distribution_plot
    File mg_distribution_plot = downstream.mg_distribution_plot
    String stat_read_count = downstream.stat_read_count
    String stat_read_length_mean = downstream.stat_read_length_mean
    String stat_read_length_median = downstream.stat_read_length_median
    String stat_read_length_n50 = downstream.stat_read_length_n50
    String stat_read_quality_mean = downstream.stat_read_quality_mean
    String stat_read_quality_median = downstream.stat_read_quality_median
    String stat_mapped_read_count = downstream.stat_mapped_read_count
    String stat_mapped_read_percent = downstream.stat_mapped_read_percent
    String stat_gap_compressed_identity_mean = downstream.stat_gap_compressed_identity_mean
    String stat_gap_compressed_identity_median = downstream.stat_gap_compressed_identity_median

    # merged, haplotagged alignments
    File merged_haplotagged_bam = downstream.merged_haplotagged_bam
    File merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index

    # mosdepth outputs
    File mosdepth_summary = upstream.mosdepth_summary
    File mosdepth_region_bed = upstream.mosdepth_region_bed
    File mosdepth_region_bed_index = upstream.mosdepth_region_bed_index
    File mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    String stat_depth_mean = upstream.stat_depth_mean
    String inferred_sex = upstream.inferred_sex

    # phasing stats
    File phase_stats = downstream.phase_stats
    File phase_blocks = downstream.phase_blocks
    File phase_haplotags = downstream.phase_haplotags
    String stat_phased_basepairs = downstream.stat_phased_basepairs
    String stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # methylation outputs and profile
    File? cpg_pileup_bed = downstream.cpg_pileup_bed
    File? cpg_pileup_bed_index = downstream.cpg_pileup_bed_index
    File? hmcpg_pileup_bed = downstream.hmcpg_pileup_bed
    File? hmcpg_pileup_bed_index = downstream.hmcpg_pileup_bed_index
    String stat_cpg_hap1_count = downstream.stat_cpg_hap1_count
    String stat_cpg_hap2_count = downstream.stat_cpg_hap2_count
    String stat_cpg_combined_count = downstream.stat_cpg_combined_count
    File? methbat_profile = downstream.methbat_profile
    String stat_methbat_methylated_count = downstream.stat_methbat_methylated_count
    String stat_methbat_unmethylated_count = downstream.stat_methbat_unmethylated_count
    String stat_methbat_asm_count = downstream.stat_methbat_asm_count

    # sv outputs
    File phased_sv_vcf = downstream.phased_sv_vcf
    File phased_sv_vcf_index = downstream.phased_sv_vcf_index
    File sv_supporting_reads = select_first([
      upstream.sv_supporting_reads
    ])
    File sv_copynum_bedgraph = select_first([
      upstream.sv_copynum_bedgraph
    ])
    File sv_depth_bw = select_first([
      upstream.sv_depth_bw
    ])
    File sv_gc_bias_corrected_depth_bw = select_first([
      upstream.sv_gc_bias_corrected_depth_bw
    ])
    File sv_copynum_summary = select_first([
      upstream.sv_copynum_summary
    ])

    # sv stats
    String stat_sv_DUP_count = downstream.stat_sv_DUP_count
    String stat_sv_DEL_count = downstream.stat_sv_DEL_count
    String stat_sv_INS_count = downstream.stat_sv_INS_count
    String stat_sv_INV_count = downstream.stat_sv_INV_count
    String stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    String stat_sv_BND_count = downstream.stat_sv_BND_count
    File sv_stats_plot = downstream.sv_stats_plot

    # small variant outputs
    File phased_small_variant_vcf = downstream.phased_small_variant_vcf
    File phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    File? small_variant_gvcf = upstream.small_variant_gvcf
    File? small_variant_gvcf_index = upstream.small_variant_gvcf_index

    # small variant stats
    File small_variant_stats = downstream.small_variant_stats
    File bcftools_roh_out = downstream.bcftools_roh_out
    File bcftools_roh_bed = downstream.bcftools_roh_bed
    String stat_small_variant_SNV_count = downstream.stat_SNV_count
    String stat_small_variant_INDEL_count = downstream.stat_INDEL_count
    String stat_small_variant_TSTV_ratio = downstream.stat_TSTV_ratio
    String stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    File snv_distribution_plot = downstream.snv_distribution_plot
    File indel_distribution_plot = downstream.indel_distribution_plot

    # trgt outputs
    File phased_trgt_vcf = downstream.trgt_vcf
    File phased_trgt_vcf_index = downstream.trgt_vcf_index
    File trgt_spanning_reads = downstream.trgt_spanning_reads
    File trgt_spanning_reads_index = downstream.trgt_spanning_reads_index
    File trgt_coverage_dropouts = downstream.trgt_coverage_dropouts
    String stat_trgt_genotyped_count = downstream.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count = downstream.stat_trgt_uncalled_count

    # paraphase outputs
    File? paraphase_summary = upstream.paraphase_output_json
    File? paraphase_realigned_bam = upstream.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    File? paraphase_vcfs = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    File mitorsaw_vcf = upstream.mitorsaw_vcf
    File mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    File mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # kivvi kiv2 outputs
    File? kivvi_kiv2_vcf = upstream.kivvi_kiv2_vcf
    File? kivvi_kiv2_vcf_index = upstream.kivvi_kiv2_vcf_index
    File? kivvi_kiv2_json = upstream.kivvi_kiv2_json
    File? kivvi_kiv2_realigned_bam = upstream.kivvi_kiv2_realigned_bam
    File? kivvi_kiv2_realigned_bam_index = upstream.kivvi_kiv2_realigned_bam_index
    File? kivvi_kiv2_allele_plot = upstream.kivvi_kiv2_allele_plot

    # kivvi d4z4 outputs
    File? kivvi_d4z4_vcf = upstream.kivvi_d4z4_vcf
    File? kivvi_d4z4_vcf_index = upstream.kivvi_d4z4_vcf_index
    File? kivvi_d4z4_json = upstream.kivvi_d4z4_json
    File? kivvi_d4z4_realigned_bam = upstream.kivvi_d4z4_realigned_bam
    File? kivvi_d4z4_realigned_bam_index = upstream.kivvi_d4z4_realigned_bam_index
    File? kivvi_d4z4_allele_plot = upstream.kivvi_d4z4_allele_plot

    # PGx outputs
    File? pbstarphase_summary = downstream.pbstarphase_json
    File? pbstarphase_tsv = downstream.pbstarphase_tsv

    # qc messages
    Array[String] msg = select_all(non_empty_msg)

    # workflow metadata
    String workflow_name = "humanwgs_singleton"
    String workflow_version = "v4.0.0-rc2" + if defined(debug_version)
      then "~{"-" + debug_version}"
      else ""
  }
}

