version 1.0

import "downstream/downstream.wdl" as Downstream
import "humanwgs_structs.wdl"
import "joint/joint.wdl" as Joint
import "process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "unpack_container_manifest/unpack_container_manifest.wdl" as UnpackContainerManifest
import "upstream/upstream.wdl" as Upstream
import "wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "wdl-common/wdl/workflows/filter_messages/filter_messages.wdl" as FilterMessages

workflow humanwgs_family {
  meta {
    version: "4.0.0"
    apiVersion: "2.1.0"
    name: "HiFi Human WGS Variant Pipeline (family)"
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
    description: "PacBio HiFi human whole genome sequencing pipeline, with joint calling for related samples."
    category: "Entrypoint"
    outputs: {
      sample_ids: {
        name: "Sample IDs",
        description: "Sample IDs"
      },
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
        name: "Read length plots",
        description: "Distribution of read lengths",
        fileTypeId: "PacBio.FileTypes.png"
      },
      read_quality_plot: {
        name: "Read quality plots",
        description: "Distribution of read qualities",
        fileTypeId: "PacBio.FileTypes.png"
      },
      mapq_distribution_plot: {
        name: "Mapping quality plots",
        description: "Distribution of mapping quality per alignment",
        fileTypeId: "PacBio.FileTypes.png"
      },
      mg_distribution_plot: {
        name: "Gap-compressed identity plots",
        description: "Distribution of gap-compressed identity per alignment",
        fileTypeId: "PacBio.FileTypes.png"
      },
      stat_read_count: {
        name: "Read counts",
        description: "Number of reads"
      },
      stat_read_length_mean: {
        name: "Mean read lengths",
        description: "Mean read length"
      },
      stat_read_length_median: {
        name: "Median read lengths",
        description: "Median read length"
      },
      stat_read_length_n50: {
        name: "Read length N50s",
        description: "Read length N50"
      },
      stat_read_quality_mean: {
        name: "Mean read qualities",
        description: "Mean read quality"
      },
      stat_read_quality_median: {
        name: "Median read qualities",
        description: "Median read quality"
      },
      stat_mapped_read_count: {
        name: "Mapped read counts",
        description: "Number of reads mapped to reference"
      },
      stat_mapped_read_percent: {
        name: "Mapped read percentages",
        description: "Percent of reads mapped to reference"
      },
      stat_gap_compressed_identity_mean: {
        name: "Mean gap-compressed identities",
        description: "Mean gap-compressed identity"
      },
      stat_gap_compressed_identity_median: {
        name: "Median gap-compressed identities",
        description: "Median gap-compressed identity"
      },
      merged_haplotagged_bam: {
        name: "Haplotagged BAMs",
        description: "Merged, haplotagged alignments",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      merged_haplotagged_bam_index: {
        name: "Haplotagged BAM indices",
        description: "Index for merged, haplotagged alignments",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      mosdepth_summary: {
        name: "Depth summaries",
        description: "Summary of aligned read depth",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      mosdepth_region_bed: {
        name: "Regional depth BEDs",
        description: "Median aligned read depth by 500bp windows",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      mosdepth_region_bed_index: {
        name: "Regional depth BED indices",
        description: "Index for median aligned read depth by 500bp windows",
        fileTypeId: "PacBio.Index.Tabix"
      },
      mosdepth_depth_distribution_plot: {
        name: "Depth distribution plots",
        description: "Distribution of aligned read depth",
        fileTypeId: "PacBio.FileTypes.png"
      },
      stat_depth_mean: {
        name: "Mean depths",
        description: "Mean depth"
      },
      inferred_sex: {
        name: "Inferred sexes",
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
        name: "Phase block NG50s",
        description: "Phase block NG50"
      },
      cpg_pileup_bed: {
        name: "5mCpG pileup BEDs",
        description: "5mCpG pileup BED",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      cpg_pileup_bed_index: {
        name: "5mCpG pileup BED indices",
        description: "Index for 5mCpG pileup BED",
        fileTypeId: "PacBio.Index.Tabix"
      },
      hmcpg_pileup_bed: {
        name: "5hmC pileup BEDs",
        description: "5hmC pileup BED",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      hmcpg_pileup_bed_index: {
        name: "5hmC pileup BED indices",
        description: "Index for 5hmC pileup BED",
        fileTypeId: "PacBio.Index.Tabix"
      },
      stat_cpg_hap1_count: {
        name: "5mCpG counts (hap1)",
        description: "Number of scored reference 5mCpGs in haplotype 1"
      },
      stat_cpg_hap2_count: {
        name: "5mCpG counts (hap2)",
        description: "Number of scored reference 5mCpGs in haplotype 2"
      },
      stat_cpg_combined_count: {
        name: "5mCpG counts (combined)",
        description: "Number of scored reference 5mCpGs combined"
      },
      methbat_profile: {
        name: "MethBat profiles",
        description: "MethBat 5mCpG profile",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      stat_methbat_methylated_count: {
        name: "Methylated region counts",
        description: "Number of profiled regions labeled as methylated"
      },
      stat_methbat_unmethylated_count: {
        name: "Unmethylated region counts",
        description: "Number of profiled regions labeled as unmethylated"
      },
      stat_methbat_asm_count: {
        name: "ASM region counts",
        description: "Number of profiled regions labeled as having allele-specific methylation"
      },
      phased_sv_vcf: {
        name: "Structural variant VCFs",
        description: "Phased structural variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      phased_sv_vcf_index: {
        name: "Structural variant VCF indices",
        description: "Index for phased structural variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      sv_supporting_reads: {
        name: "SV supporting reads",
        description: "Supporting reads for structural variants",
        fileTypeId: "PacBio.FileTypes.json"
      },
      sv_copynum_bedgraph: {
        name: "CNV copy number BEDGraphs",
        description: "CNV copy number BEDGraph",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      sv_depth_bw: {
        name: "CNV depth BigWigs",
        description: "CNV depth BigWig",
        fileTypeId: "PacBio.FileTypes.bigwig"
      },
      sv_gc_bias_corrected_depth_bw: {
        name: "CNV GC-corrected depth BigWigs",
        description: "CNV GC-bias corrected depth BigWig",
        fileTypeId: "PacBio.FileTypes.bigwig"
      },
      sv_copynum_summary: {
        name: "CNV copy number summaries",
        description: "CNV copy number summary JSON",
        fileTypeId: "PacBio.FileTypes.json"
      },
      stat_sv_DUP_count: {
        name: "Duplication counts",
        description: "Number of DUP structural variants"
      },
      stat_sv_DEL_count: {
        name: "Deletion counts",
        description: "Number of DEL structural variants"
      },
      stat_sv_INS_count: {
        name: "Insertion counts",
        description: "Number of INS structural variants"
      },
      stat_sv_INV_count: {
        name: "Inversion counts",
        description: "Number of INV structural variants"
      },
      stat_sv_SWAP_count: {
        name: "Sequence swap counts",
        description: "Number of structural variant sequence swap events"
      },
      stat_sv_BND_count: {
        name: "Breakend counts",
        description: "Number of BND structural variants"
      },
      sv_stats_plot: {
        name: "SV size distribution plots",
        description: "Structural variant size distribution plot",
        fileTypeId: "PacBio.FileTypes.png"
      },
      phased_small_variant_vcf: {
        name: "Small variant VCFs",
        description: "Phased small variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      phased_small_variant_vcf_index: {
        name: "Small variant VCF indices",
        description: "Index for phased small variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      small_variant_gvcf: {
        name: "Small variant GVCFs",
        description: "Small variant GVCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      small_variant_gvcf_index: {
        name: "Small variant GVCF indices",
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
        name: "Regions of homozygosity BEDs",
        description: "Regions of homozygosity BED",
        fileTypeId: "PacBio.FileTypes.bed"
      },
      stat_small_variant_SNV_count: {
        name: "SNV counts",
        description: "Number of SNVs"
      },
      stat_small_variant_INDEL_count: {
        name: "INDEL counts",
        description: "Number of INDELs"
      },
      stat_small_variant_TSTV_ratio: {
        name: "Ts/Tv ratios",
        description: "Ts/Tv ratio"
      },
      stat_small_variant_HETHOM_ratio: {
        name: "Het/Hom ratios",
        description: "Het/Hom ratio for SNVs"
      },
      snv_distribution_plot: {
        name: "SNV distribution plots",
        description: "Distribution of SNVs by REF, ALT",
        fileTypeId: "PacBio.FileTypes.png"
      },
      indel_distribution_plot: {
        name: "INDEL distribution plots",
        description: "Distribution of indels by size",
        fileTypeId: "PacBio.FileTypes.png"
      },
      phased_trgt_vcf: {
        name: "Tandem repeat VCFs",
        description: "Phased TRGT VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      phased_trgt_vcf_index: {
        name: "Tandem repeat VCF indices",
        description: "Index for phased TRGT VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      trgt_spanning_reads: {
        name: "Tandem repeat spanning reads",
        description: "Aligned TRGT spanning reads",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      trgt_spanning_reads_index: {
        name: "Spanning reads indices",
        description: "Index for aligned TRGT spanning reads",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      trgt_coverage_dropouts: {
        name: "Coverage dropout regions",
        description: "TRGT regions with coverage dropouts",
        fileTypeId: "PacBio.FileTypes.txt"
      },
      stat_trgt_genotyped_count: {
        name: "Genotyped repeat counts",
        description: "Number of sites genotyped by TRGT"
      },
      stat_trgt_uncalled_count: {
        name: "Ungenotyped repeat counts",
        description: "Number of sites ungenotyped by TRGT"
      },
      paraphase_summary: {
        name: "Paraphase summaries",
        description: "Paraphase summary",
        fileTypeId: "PacBio.FileTypes.json"
      },
      paraphase_realigned_bam: {
        name: "Paraphase-realigned BAMs",
        description: "BAM file of reads realigned by Paraphase",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      paraphase_realigned_bam_index: {
        name: "Paraphase-realigned BAM indices",
        description: "Index for BAM file of reads realigned by Paraphase",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      paraphase_vcfs: {
        name: "Paraphase VCFs",
        description: "Paraphase VCFs",
        fileTypeId: "PacBio.FileTypes.tgz"
      },
      mitorsaw_vcf: {
        name: "Mitochondrial variant VCFs",
        description: "Mitochondrial variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      mitorsaw_vcf_index: {
        name: "Mitochondrial variant VCF indices",
        description: "Index for mitochondrial variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      mitorsaw_hap_stats: {
        name: "Mitochondrial haplotype statistics",
        description: "Mitochondrial haplotype statistics",
        fileTypeId: "PacBio.FileTypes.json"
      },
      kivvi_kiv2_vcf: {
        name: "KIV2 variant VCFs",
        description: "KIV2 repeat variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      kivvi_kiv2_vcf_index: {
        name: "KIV2 variant VCF indices",
        description: "Index for KIV2 repeat variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      kivvi_kiv2_json: {
        name: "KIV2 genotype JSONs",
        description: "KIV2 repeat genotype JSON",
        fileTypeId: "PacBio.FileTypes.json"
      },
      kivvi_kiv2_realigned_bam: {
        name: "KIV2-realigned BAMs",
        description: "KIV2-realigned BAM",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      kivvi_kiv2_realigned_bam_index: {
        name: "KIV2-realigned BAM indices",
        description: "Index for KIV2-realigned BAM",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      kivvi_kiv2_allele_plot: {
        name: "KIV2 allele plots",
        description: "KIV2 assembled allele plot",
        fileTypeId: "PacBio.FileTypes.svg"
      },
      kivvi_d4z4_vcf: {
        name: "D4Z4 variant VCFs",
        description: "D4Z4 repeat variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      kivvi_d4z4_vcf_index: {
        name: "D4Z4 variant VCF indices",
        description: "Index for D4Z4 repeat variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      kivvi_d4z4_json: {
        name: "D4Z4 genotype JSONs",
        description: "D4Z4 repeat genotype JSON",
        fileTypeId: "PacBio.FileTypes.json"
      },
      kivvi_d4z4_realigned_bam: {
        name: "D4Z4-realigned BAMs",
        description: "D4Z4-realigned BAM",
        fileTypeId: "PacBio.AlignmentFile.ConsensusAlignmentBamFile"
      },
      kivvi_d4z4_realigned_bam_index: {
        name: "D4Z4-realigned BAM indices",
        description: "Index for D4Z4-realigned BAM",
        fileTypeId: "PacBio.Index.BamIndex"
      },
      kivvi_d4z4_allele_plot: {
        name: "D4Z4 allele plots",
        description: "D4Z4 assembled allele plot",
        fileTypeId: "PacBio.FileTypes.svg"
      },
      pbstarphase_summary: {
        name: "StarPhase summaries",
        description: "StarPhase summary",
        fileTypeId: "PacBio.FileTypes.json"
      },
      pbstarphase_tsv: {
        name: "StarPhase PharmCAT TSVs",
        description: "StarPhase summary in TSV format for PharmCAT",
        fileTypeId: "PacBio.FileTypes.tsv"
      },
      joint_small_variants_vcf: {
        name: "Joint-called small variant VCF",
        description: "Joint-called small variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      joint_small_variants_vcf_index: {
        name: "Joint-called small variant VCF index",
        description: "Index for joint-called small variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      joint_sv_vcf: {
        name: "Joint-called structural variant VCF",
        description: "Joint-called structural variant VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      joint_sv_vcf_index: {
        name: "Joint-called structural variant VCF index",
        description: "Index for joint-called structural variant VCF",
        fileTypeId: "PacBio.Index.Tabix"
      },
      joint_trgt_vcf: {
        name: "Joint-called tandem repeat VCF",
        description: "Joint-called TRGT VCF",
        fileTypeId: "PacBio.FileTypes.vcf"
      },
      joint_trgt_vcf_index: {
        name: "Joint-called tandem repeat VCF index",
        description: "Index for joint-called TRGT VCF",
        fileTypeId: "PacBio.Index.Tabix"
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
    family: {
      name: "Family",
      description: "Family struct describing samples, relationships, and unaligned BAM paths",
      advanced: false
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
    glnexus_mem_gb: {
      name: "GLnexus memory override",
      description: "Override GLnexus memory request (GB)",
      group: "Resource"
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
  }

  input {
    Family family
    String ref_name = "GRCh38_GIABv3"
    File? trgt_tandem_repeat_bed_override
    File? methbat_region_tsv_override
    Boolean use_alignment_chunking = true
    Int glnexus_mem_gb = 60
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

  Boolean single_sample = length(family.samples) == 1

  scatter (sample in family.samples) {
    String sample_id = sample.sample_id

    call Upstream.upstream { input:
      sample_id = sample.sample_id,
      hifi_reads = sample.hifi_reads,
      fail_reads = sample.fail_reads,
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
      single_sample = single_sample,
      use_gpu = use_gpu,
      use_parabricks_deepvariant = use_parabricks_deepvariant,
      default_runtime_attributes = default_runtime_attributes
    }
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
      ref_name = ref_name,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
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
      ref_name = ref_name,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      sawfish_expected_bed_male = sawfish_expected_bed_male,
      sawfish_expected_bed_female = sawfish_expected_bed_female,
      methbat_region_tsv = methbat_region_tsv,
      run_starphase = run_starphase,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  if (!single_sample) {
    call Bcftools.bcftools_merge as merge_small_variant_vcfs { input:
      vcfs = downstream.phased_small_variant_vcf,
      vcf_indices = downstream.phased_small_variant_vcf_index,
      out_prefix = "~{family.family_id}.joint.~{ref_name}.small_variants.phased",
      runtime_attributes = default_runtime_attributes
    }

    call Bcftools.bcftools_merge as merge_sv_vcfs { input:
      vcfs = downstream.phased_sv_vcf,
      vcf_indices = downstream.phased_sv_vcf_index,
      out_prefix = "~{family.family_id}.joint.~{ref_name}.structural_variants.phased",
      runtime_attributes = default_runtime_attributes
    }

    call Trgt.trgt_merge { input:
      vcfs = downstream.trgt_vcf,
      vcf_indices = downstream.trgt_vcf_index,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      out_prefix = "~{family.family_id}.merged.~{ref_name}.trgt",
      runtime_attributes = default_runtime_attributes
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
      downstream.stat_cpg_combined_count
    ]),
    flatten([
      [
        "cpg_hap1_count"
      ],
      downstream.stat_cpg_hap1_count
    ]),
    flatten([
      [
        "cpg_hap2_count"
      ],
      downstream.stat_cpg_hap2_count
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

  call FilterMessages.filter_messages { input:
    message_arrays = [
      process_trgt_catalog.msg,
      flatten(upstream.msg),
      flatten(downstream.msg)
    ]
  }

  call Utilities.consolidate_stats { input:
    out_prefix = family.family_id,
    stats = stats,
    msg_array = filter_messages.messages,
    runtime_attributes = default_runtime_attributes
  }

  output {
    # to maintain order of samples
    Array[String] sample_ids = sample_id
    File stats_file = consolidate_stats.stats_tsv
    File msg_file = consolidate_stats.messages

    # bam stats
    Array[File] read_length_plot = downstream.read_length_plot
    Array[File] read_quality_plot = downstream.read_quality_plot
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
    Array[File?] cpg_pileup_bed = downstream.cpg_pileup_bed
    Array[File?] cpg_pileup_bed_index = downstream.cpg_pileup_bed_index
    Array[File?] hmcpg_pileup_bed = downstream.hmcpg_pileup_bed
    Array[File?] hmcpg_pileup_bed_index = downstream.hmcpg_pileup_bed_index
    Array[String] stat_cpg_hap1_count = downstream.stat_cpg_hap1_count
    Array[String] stat_cpg_hap2_count = downstream.stat_cpg_hap2_count
    Array[String] stat_cpg_combined_count = downstream.stat_cpg_combined_count
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
    Array[File] sv_stats_plot = downstream.sv_stats_plot

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

    # per sample kivvi kiv2 outputs
    Array[File?] kivvi_kiv2_vcf = upstream.kivvi_kiv2_vcf
    Array[File?] kivvi_kiv2_vcf_index = upstream.kivvi_kiv2_vcf_index
    Array[File?] kivvi_kiv2_json = upstream.kivvi_kiv2_json
    Array[File?] kivvi_kiv2_realigned_bam = upstream.kivvi_kiv2_realigned_bam
    Array[File?] kivvi_kiv2_realigned_bam_index = upstream.kivvi_kiv2_realigned_bam_index
    Array[File?] kivvi_kiv2_allele_plot = upstream.kivvi_kiv2_allele_plot

    # per sample kivvi d4z4 outputs
    Array[File?] kivvi_d4z4_vcf = upstream.kivvi_d4z4_vcf
    Array[File?] kivvi_d4z4_vcf_index = upstream.kivvi_d4z4_vcf_index
    Array[File?] kivvi_d4z4_json = upstream.kivvi_d4z4_json
    Array[File?] kivvi_d4z4_realigned_bam = upstream.kivvi_d4z4_realigned_bam
    Array[File?] kivvi_d4z4_realigned_bam_index = upstream.kivvi_d4z4_realigned_bam_index
    Array[File?] kivvi_d4z4_allele_plot = upstream.kivvi_d4z4_allele_plot

    # PGx outputs
    Array[File?] pbstarphase_summary = downstream.pbstarphase_json
    Array[File?] pbstarphase_tsv = downstream.pbstarphase_tsv

    # joint call outputs
    File? joint_small_variants_vcf = merge_small_variant_vcfs.merged_vcf
    File? joint_small_variants_vcf_index = merge_small_variant_vcfs.merged_vcf_index
    File? joint_sv_vcf = merge_sv_vcfs.merged_vcf
    File? joint_sv_vcf_index = merge_sv_vcfs.merged_vcf_index
    File? joint_trgt_vcf = trgt_merge.merged_vcf
    File? joint_trgt_vcf_index = trgt_merge.merged_vcf_index

    # qc messages
    Array[String] msg = filter_messages.messages

    # workflow metadata
    String workflow_name = "humanwgs_family"
    String workflow_version = "4.0.0" + if defined(debug_version)
      then "~{"-" + debug_version}"
      else ""
  }
}

