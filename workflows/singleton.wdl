version 1.0

import "downstream/downstream.wdl" as Downstream
import "process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "unpack_container_manifest/unpack_container_manifest.wdl" as UnpackContainerManifest
import "upstream/upstream.wdl" as Upstream
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "wdl-common/wdl/workflows/filter_messages/filter_messages.wdl" as FilterMessages

workflow humanwgs_singleton {
  meta {
    version: "4.0.0"
    apiVersion: "2.1.0"
    name: "HiFi Human WGS Variant Pipeline (singleton)"
    isAutoAnalysis: true
    supportedBackends: [
      "HPC",
      "AWS-HealthOmics"
    ]
    tags: "ccs analysis auto-analysis"
    containers: [
      "quay.io/pacbio/glnexus@sha256:ce6fecf59dddc6089a8100b31c29c1e6ed50a0cf123da9f2bc589ee4b0c69c8e",
      "quay.io/pacbio/hiphase@sha256:41ebe22b55c66e2e78da2013f7fffaecc02a8b4e980400c3ea8d03c87330522e",
      "quay.io/pacbio/kivvi@sha256:9e4a390821ea999af9b3faa483c5f7dce58041c439a35e4af3676f25edecd204",
      "quay.io/pacbio/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb",
      "quay.io/pacbio/mitorsaw@sha256:9b610f343ae018b21c86977259d91dba587e1f22d7a73e115937e428d3c7adeb",
      "quay.io/pacbio/mosdepth@sha256:f4edf52e6a31eb18f1755e8cb90224fec6aa88e4a4353a673b16a89f59cbe822",
      "quay.io/pacbio/paraphase@sha256:665ac4fefcef92e0023395eae9e0ca3f6e65d1228c49afba3ad3f8e7b3f3d1cb",
      "quay.io/pacbio/pbjam@sha256:1b90537ed6683ac7b89002c0f59e940b4ea46c600d292904f549af1ae2760bba",
      "quay.io/pacbio/pbmm2@sha256:86e39aa67fa5d385769d5f119739e8811ce550163e3b9dfc42bd58d1fecdf3a8",
      "quay.io/pacbio/pbsamoa@sha256:cf70c89422d63c3a4e4b2ffd912160c3e7481303a211cd6134236fe6e9f9461f",
      "quay.io/pacbio/pbstarphase@sha256:86aeeb3a9663c22c135d08353de700050e1d268406c1f63832cad5f139cdd390",
      "quay.io/pacbio/pbtk@sha256:f27bafa0ae6ffff6170d45a86a5677406320c7866760391f89ebe900a74d2039",
      "quay.io/pacbio/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae",
      "quay.io/pacbio/sawfish@sha256:5cf8f02790ecb89e652885c57f103fc9597c41f75ee71be2c05510cbbdf68f59",
      "quay.io/pacbio/trgt@sha256:648aee4a2c9d7371a48e454a7143861a242b853d81ff5453924cc0095d207824",
      "google/deepvariant:1.10.0",
      "google/deepvariant:1.10.0-gpu",
      "nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1",
      "quay.io/pacbio/workflow-data-container-hifi-human-wgs-wdl-grch38@sha256:5e3f23e44c1c09762838e81677f7d136367be2a4063efc608129841cef745b42",
      "quay.io/pacbio/workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3@sha256:f3799c1eff1b816a95ea06ce567d3324d06e945428567077d96eca07c76ad0aa"
    ]
    description: "PacBio HiFi human whole genome sequencing pipeline for individual samples."
    category: "Entrypoint"
    outputs: {
      stats_file: {
        name: "Summary statistics",
        description: "Table of summary statistics",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      msg_file: {
        name: "Messages file",
        description: "File containing messages from the workflow",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      read_length_plot: {
        name: "Read length plot",
        description: "Distribution of read lengths",
        fileTypeId: "PacBio.FileTypes.png"
      },
      read_quality_plot: {
        name: "Read quality plot",
        description: "Distribution of read qualities",
        fileTypeId: "PacBio.FileTypes.png"
      },
      mapq_distribution_plot: {
        name: "Mapping quality plot",
        description: "Distribution of mapping quality per alignment",
        fileTypeId: "PacBio.FileTypes.png"
      },
      mg_distribution_plot: {
        name: "Gap-compressed identity plot",
        description: "Distribution of gap-compressed identity per alignment",
        fileTypeId: "PacBio.FileTypes.png"
      },
      stat_read_count: {
        name: "Read count",
        description: "Number of reads"
      },
      stat_read_length_mean: {
        name: "Mean read length",
        description: "Mean read length"
      },
      stat_read_length_median: {
        name: "Median read length",
        description: "Median read length"
      },
      stat_read_length_n50: {
        name: "Read length N50",
        description: "Read length N50"
      },
      stat_read_quality_mean: {
        name: "Mean read quality",
        description: "Mean read quality"
      },
      stat_read_quality_median: {
        name: "Median read quality",
        description: "Median read quality"
      },
      stat_mapped_read_count: {
        name: "Mapped read count",
        description: "Number of reads mapped to reference"
      },
      stat_mapped_read_percent: {
        name: "Mapped read percent",
        description: "Percent of reads mapped to reference"
      },
      stat_gap_compressed_identity_mean: {
        name: "Mean gap-compressed identity",
        description: "Mean gap-compressed identity"
      },
      stat_gap_compressed_identity_median: {
        name: "Median gap-compressed identity",
        description: "Median gap-compressed identity"
      },
      merged_haplotagged_bam: {
        name: "Haplotagged BAM",
        description: "Merged, haplotagged alignments",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      merged_haplotagged_bam_index: {
        name: "Haplotagged BAM index",
        description: "Index for merged, haplotagged alignments",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      mosdepth_summary: {
        name: "Depth summary",
        description: "Summary of aligned read depth",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      mosdepth_region_bed: {
        name: "Regional depth BED",
        description: "Median aligned read depth by 500bp windows",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      mosdepth_region_bed_index: {
        name: "Regional depth BED index",
        description: "Index for median aligned read depth by 500bp windows",
        fileTypeId: "PacBio.Index.Tabix"
      },
      mosdepth_depth_distribution_plot: {
        name: "Depth distribution plot",
        description: "Distribution of aligned read depth",
        fileTypeId: "PacBio.FileTypes.png"
      },
      stat_depth_mean: {
        name: "Mean depth",
        description: "Mean depth"
      },
      inferred_sex: {
        name: "Inferred sex",
        description: "Inferred sex"
      },
      phase_stats: {
        name: "Phasing statistics",
        description: "Phasing statistics",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      phase_blocks: {
        name: "Phase blocks",
        description: "Phase blocks",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      phase_haplotags: {
        name: "Read phase assignments",
        description: "Per-read phase assignment",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      stat_phased_basepairs: {
        name: "Phased basepairs",
        description: "Number of basepairs within phase blocks"
      },
      stat_phase_block_ng50: {
        name: "Phase block NG50",
        description: "Phase block NG50"
      },
      cpg_pileup_bed: {
        name: "5mCpG pileup BED",
        description: "5mCpG pileup BED",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      cpg_pileup_bed_index: {
        name: "5mCpG pileup BED index",
        description: "Index for 5mCpG pileup BED",
        fileTypeId: "PacBio.Index.Tabix"
      },
      hmcpg_pileup_bed: {
        name: "5hmC pileup BED",
        description: "5hmC pileup BED",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      hmcpg_pileup_bed_index: {
        name: "5hmC pileup BED index",
        description: "Index for 5hmC pileup BED",
        fileTypeId: "PacBio.Index.Tabix"
      },
      stat_cpg_hap1_count: {
        name: "5mCpG count (hap1)",
        description: "Number of scored reference 5mCpGs in haplotype 1"
      },
      stat_cpg_hap2_count: {
        name: "5mCpG count (hap2)",
        description: "Number of scored reference 5mCpGs in haplotype 2"
      },
      stat_cpg_combined_count: {
        name: "5mCpG count (combined)",
        description: "Number of scored reference 5mCpGs combined"
      },
      methbat_profile: {
        name: "MethBat profile",
        description: "MethBat 5mCpG profile",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      stat_methbat_methylated_count: {
        name: "Methylated region count",
        description: "Number of profiled regions labeled as methylated"
      },
      stat_methbat_unmethylated_count: {
        name: "Unmethylated region count",
        description: "Number of profiled regions labeled as unmethylated"
      },
      stat_methbat_asm_count: {
        name: "ASM region count",
        description: "Number of profiled regions labeled as having allele-specific methylation"
      },
      phased_sv_vcf: {
        name: "Structural variant VCF",
        description: "Phased structural variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      phased_sv_vcf_index: {
        name: "Structural variant VCF index",
        description: "Index for phased structural variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      sv_supporting_reads: {
        name: "SV supporting reads",
        description: "Supporting reads for structural variants",
        fileTypeId: "PacBio.FileTypes.json"
      },
      sv_copynum_bedgraph: {
        name: "CNV copy number BEDGraph",
        description: "CNV copy number BEDGraph",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      sv_depth_bw: {
        name: "CNV depth BigWig",
        description: "CNV depth BigWig",
        fileTypeId: "PacBio.FileTypes.bigwig"
      },
      sv_gc_bias_corrected_depth_bw: {
        name: "CNV GC-corrected depth BigWig",
        description: "CNV GC-bias corrected depth BigWig",
        fileTypeId: "PacBio.FileTypes.bigwig"
      },
      sv_copynum_summary: {
        name: "CNV copy number summary",
        description: "CNV copy number summary JSON",
        fileTypeId: "PacBio.FileTypes.json"
      },
      stat_sv_DUP_count: {
        name: "Duplication count",
        description: "Number of DUP structural variants"
      },
      stat_sv_DEL_count: {
        name: "Deletion count",
        description: "Number of DEL structural variants"
      },
      stat_sv_INS_count: {
        name: "Insertion count",
        description: "Number of INS structural variants"
      },
      stat_sv_INV_count: {
        name: "Inversion count",
        description: "Number of INV structural variants"
      },
      stat_sv_SWAP_count: {
        name: "Sequence swap count",
        description: "Number of structural variant sequence swap events"
      },
      stat_sv_BND_count: {
        name: "Breakend count",
        description: "Number of BND structural variants"
      },
      sv_stats_plot: {
        name: "SV size distribution plot",
        description: "Structural variant size distribution plot",
        fileTypeId: "PacBio.FileTypes.png"
      },
      phased_small_variant_vcf: {
        name: "Small variant VCF",
        description: "Phased small variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      phased_small_variant_vcf_index: {
        name: "Small variant VCF index",
        description: "Index for phased small variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      small_variant_gvcf: {
        name: "Small variant GVCF",
        description: "Small variant GVCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      small_variant_gvcf_index: {
        name: "Small variant GVCF index",
        description: "Index for small variant GVCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      small_variant_stats: {
        name: "Small variant statistics",
        description: "Small variant statistics",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      bcftools_roh_out: {
        name: "Regions of homozygosity",
        description: "Regions of homozygosity",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      bcftools_roh_bed: {
        name: "Regions of homozygosity BED",
        description: "Regions of homozygosity BED",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      stat_small_variant_SNV_count: {
        name: "SNV count",
        description: "Number of SNVs"
      },
      stat_small_variant_INDEL_count: {
        name: "INDEL count",
        description: "Number of INDELs"
      },
      stat_small_variant_TSTV_ratio: {
        name: "Ts/Tv ratio",
        description: "Ts/Tv ratio"
      },
      stat_small_variant_HETHOM_ratio: {
        name: "Het/Hom ratio",
        description: "Het/Hom ratio for SNVs"
      },
      snv_distribution_plot: {
        name: "SNV distribution plot",
        description: "Distribution of SNVs by REF, ALT",
        fileTypeId: "PacBio.FileTypes.png"
      },
      indel_distribution_plot: {
        name: "INDEL distribution plot",
        description: "Distribution of indels by size",
        fileTypeId: "PacBio.FileTypes.png"
      },
      phased_trgt_vcf: {
        name: "Tandem repeat VCF",
        description: "Phased TRGT VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      phased_trgt_vcf_index: {
        name: "Tandem repeat VCF index",
        description: "Index for phased TRGT VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      trgt_spanning_reads: {
        name: "Tandem repeat spanning reads",
        description: "Aligned TRGT spanning reads",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      trgt_spanning_reads_index: {
        name: "Spanning reads index",
        description: "Index for aligned TRGT spanning reads",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      trgt_coverage_dropouts: {
        name: "Coverage dropout regions",
        description: "TRGT regions with coverage dropouts",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      stat_trgt_genotyped_count: {
        name: "Genotyped repeat count",
        description: "Number of sites genotyped by TRGT"
      },
      stat_trgt_uncalled_count: {
        name: "Ungenotyped repeat count",
        description: "Number of sites ungenotyped by TRGT"
      },
      paraphase_summary: {
        name: "Paraphase summary",
        description: "Paraphase summary",
        fileTypeId: "PacBio.FileTypes.json"
      },
      paraphase_realigned_bam: {
        name: "Paraphase-realigned BAM",
        description: "BAM file of reads realigned by Paraphase",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      paraphase_realigned_bam_index: {
        name: "Paraphase-realigned BAM index",
        description: "Index for BAM file of reads realigned by Paraphase",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      paraphase_vcfs: {
        name: "Paraphase VCFs",
        description: "Paraphase VCFs",
        fileTypeId: "PacBio.FileTypes.tgz"
      },
      mitorsaw_vcf: {
        name: "Mitochondrial variant VCF",
        description: "Mitochondrial variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      mitorsaw_vcf_index: {
        name: "Mitochondrial variant VCF index",
        description: "Index for mitochondrial variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      mitorsaw_hap_stats: {
        name: "Mitochondrial haplotype statistics",
        description: "Mitochondrial haplotype statistics",
        fileTypeId: "PacBio.FileTypes.json"
      },
      kivvi_kiv2_vcf: {
        name: "KIV2 variant VCF",
        description: "KIV2 repeat variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      kivvi_kiv2_vcf_index: {
        name: "KIV2 variant VCF index",
        description: "Index for KIV2 repeat variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      kivvi_kiv2_json: {
        name: "KIV2 genotype JSON",
        description: "KIV2 repeat genotype JSON",
        fileTypeId: "PacBio.FileTypes.json"
      },
      kivvi_kiv2_realigned_bam: {
        name: "KIV2-realigned BAM",
        description: "KIV2-realigned BAM",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      kivvi_kiv2_realigned_bam_index: {
        name: "KIV2-realigned BAM index",
        description: "Index for KIV2-realigned BAM",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      kivvi_kiv2_allele_plot: {
        name: "KIV2 allele plot",
        description: "KIV2 assembled allele plot",
        fileTypeId: "PacBio.FileTypes.svg"
      },
      kivvi_d4z4_vcf: {
        name: "D4Z4 variant VCF",
        description: "D4Z4 repeat variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      kivvi_d4z4_vcf_index: {
        name: "D4Z4 variant VCF index",
        description: "Index for D4Z4 repeat variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      kivvi_d4z4_json: {
        name: "D4Z4 genotype JSON",
        description: "D4Z4 repeat genotype JSON",
        fileTypeId: "PacBio.FileTypes.json"
      },
      kivvi_d4z4_realigned_bam: {
        name: "D4Z4-realigned BAM",
        description: "D4Z4-realigned BAM",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      kivvi_d4z4_realigned_bam_index: {
        name: "D4Z4-realigned BAM index",
        description: "Index for D4Z4-realigned BAM",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      kivvi_d4z4_allele_plot: {
        name: "D4Z4 allele plot",
        description: "D4Z4 assembled allele plot",
        fileTypeId: "PacBio.FileTypes.svg"
      },
      pbstarphase_summary: {
        name: "StarPhase summary",
        description: "StarPhase summary",
        fileTypeId: "PacBio.FileTypes.json"
      },
      pbstarphase_tsv: {
        name: "StarPhase PharmCAT TSV",
        description: "StarPhase summary in TSV format for PharmCAT",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      msg: {
        name: "Workflow messages",
        description: "Messages from the workflow"
      },
      workflow_name: {
        name: "Workflow name",
        description: "Workflow name"
      },
      workflow_version: {
        name: "Workflow version",
        description: "Workflow version"
      }
    }
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID",
      description: "Unique identifier for the sample",
      advanced: false
    }
    hifi_reads: {
      name: "HiFi reads",
      description: "Array of paths to hifi_reads in unaligned BAM format",
      datasetType: "PacBio.DataSet.ConsensusReadSet",
      fileTypeId: "PacBio.ConsensusReadFile.ConsensusReadBamFile"
    }
    fail_reads: {
      name: "Fail reads",
      description: "Array of paths to fail_reads in unaligned BAM format",
      datasetType: "PacBio.DataSet.ConsensusReadSet",
      fileTypeId: "PacBio.ConsensusReadFile.FailReadBamFile",
      group: "Common"
    }
    ref_name: {
      name: "Reference genome",
      description: "Reference genome to use for this workflow run",
      choices: [
        "GRCh38",
        "GRCh38_GIABv3"
      ],
      default: "GRCh38_GIABv3",
      advanced: false,
      group: "Common"
    }
    trgt_tandem_repeat_bed_override: {
      name: "TRGT tandem repeat catalog override",
      description: "Optional BED file to override the default TRGT tandem repeat catalog",
      fileTypeId: "PacBio.FileTypes.bed"
    }
    methbat_region_tsv_override: {
      name: "MethBat region catalog override",
      description: "Optional TSV file to override the default MethBat methylation profiling regions",
      fileTypeId: "PacBio.FileTypes.tsv"
    }
    use_alignment_chunking: {
      name: "Use alignment chunking",
      description: "Whether to chunk BAM files for alignment. If false, all reads will be aligned in a single chunk."
    }
    use_gpu: {
      name: "Use GPU",
      description: "Use GPU when possible",
      group: "Resource"
    }
    use_parabricks_deepvariant: {
      name: "Use Parabricks DeepVariant",
      description: "Use Parabricks DeepVariant for small variant calling when GPU is enabled",
      group: "Resource"
    }
    backend: {
      name: "Backend",
      description: "Backend where the workflow will be executed",
      choices: [
        "GCP",
        "Azure",
        "AWS-HealthOmics",
        "HPC"
      ],
      group: "Resource"
    }
    zones: {
      name: "Zones",
      description: "Zones where compute will take place; required if backend is set to 'GCP'",
      hidden: true,
      group: "Resource"
    }
    cpuPlatform: {
      name: "CPU platform",
      description: "Optional minimum CPU platform to use for tasks on GCP",
      hidden: true,
      group: "Resource"
    }
    gpuType: {
      name: "GPU type",
      description: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend",
      hidden: true,
      group: "Resource"
    }
    container_registry: {
      name: "Container registry",
      description: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'",
      default: "quay.io/pacbio",
      hidden: true,
      group: "Resource"
    }
    preemptible: {
      name: "Preemptible",
      description: "Where possible, run tasks preemptibly",
      hidden: true,
      group: "Resource"
    }
    debug_version: {
      name: "Debug version",
      description: "Debug version for testing purposes",
      hidden: true
    }
    eid_ccs: {
      fileTypeId: "PacBio.DataSet.ConsensusReadSet",
      name: "HiFi Reads",
      fieldNames: [
        "hifi_reads",
        "fail_reads"
      ]
    }
  }

  input {
    String sample_id
    Array[File] hifi_reads
    Array[File]? fail_reads
    String ref_name = "GRCh38_GIABv3"
    File? trgt_tandem_repeat_bed_override
    File? methbat_region_tsv_override
    Boolean use_alignment_chunking = true
    Boolean use_gpu = false
    Boolean use_parabricks_deepvariant = false

    # Backend configuration
    String backend = "HPC"
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
  # here when a new reference is supported.
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

  call FilterMessages.filter_messages { input:
    message_arrays = [
      process_trgt_catalog.msg,
      upstream.msg,
      downstream.msg
    ]
  }

  call Utilities.consolidate_stats { input:
    out_prefix = sample_id,
    stats = stats,
    msg_array = filter_messages.messages,
    runtime_attributes = default_runtime_attributes
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
    Array[String] msg = filter_messages.messages

    # workflow metadata
    String workflow_name = "humanwgs_singleton"
    String workflow_version = "4.0.0" + if defined(debug_version)
      then "~{"-" + debug_version}"
      else ""
  }
}

