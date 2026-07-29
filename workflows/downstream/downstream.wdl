version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "../wdl-common/wdl/tasks/hiphase.wdl" as Hiphase
import "../wdl-common/wdl/tasks/methbat.wdl" as Methbat
import "../wdl-common/wdl/tasks/pbjam.wdl" as Pbjam
import "../wdl-common/wdl/tasks/pbstarphase.wdl" as Pbstarphase
import "../wdl-common/wdl/tasks/trgt.wdl" as Trgt

workflow downstream {
  meta {
    description: "Phases small variants, SVs, and TRGTs, haplotags alignments, calls HLA and PGx alleles."
    outputs: {
      merged_haplotagged_bam: {
        description: "Merged, haplotagged alignments"
      },
      merged_haplotagged_bam_index: {
        description: "Index for merged, haplotagged alignments"
      },
      phased_small_variant_vcf: {
        description: "Phased small variant VCF"
      },
      phased_small_variant_vcf_index: {
        description: "Index for phased small variant VCF"
      },
      phased_sv_vcf: {
        description: "Phased structural variant VCF"
      },
      phased_sv_vcf_index: {
        description: "Index for phased structural variant VCF"
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
      trgt_coverage_dropouts: {
        description: "TRGT regions with coverage dropouts"
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
      stat_SNV_count: {
        description: "Number of SNVs"
      },
      stat_INDEL_count: {
        description: "Number of INDELs"
      },
      stat_TSTV_ratio: {
        description: "Ts/Tv ratio"
      },
      stat_HETHOM_ratio: {
        description: "Het/Hom ratio for SNVs"
      },
      snv_distribution_plot: {
        description: "Distribution of SNVs by REF, ALT"
      },
      indel_distribution_plot: {
        description: "Distribution of indels by size"
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
      stat_sv_BND_count: {
        description: "Number of BND structural variants"
      },
      stat_sv_SWAP_count: {
        description: "Number of structural variant sequence swap events"
      },
      sv_stats_plot: {
        description: "Structural variant size distribution plot"
      },
      trgt_vcf: {
        description: "Phased TRGT VCF"
      },
      trgt_vcf_index: {
        description: "Index for phased TRGT VCF"
      },
      trgt_spanning_reads: {
        description: "Aligned TRGT spanning reads"
      },
      trgt_spanning_reads_index: {
        description: "Index for aligned TRGT spanning reads"
      },
      stat_trgt_genotyped_count: {
        description: "Number of sites genotyped by TRGT"
      },
      stat_trgt_uncalled_count: {
        description: "Number of sites ungenotyped by TRGT"
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
      pbstarphase_json: {
        description: "StarPhase summary"
      },
      pbstarphase_tsv: {
        description: "StarPhase summary in TSV format for PharmCAT"
      },
      msg: {
        description: "Messages from the workflow"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    sex: {
      description: "Sample sex",
      choices: [
        "MALE",
        "FEMALE"
      ]
    }
    aligned_hifi_reads: {
      description: "Aligned hifi_reads BAM"
    }
    aligned_hifi_reads_index: {
      description: "Aligned hifi_reads BAM index"
    }
    aligned_fail_reads: {
      description: "Aligned fail_reads BAM"
    }
    aligned_fail_reads_index: {
      description: "Aligned fail_reads BAM index"
    }
    trgt_catalog: {
      description: "TRGT tandem repeat catalog"
    }
    small_variant_vcf: {
      description: "Small variant VCF"
    }
    small_variant_vcf_index: {
      description: "Small variant VCF index"
    }
    sv_vcf: {
      description: "Structural variant VCF"
    }
    sv_vcf_index: {
      description: "Structural variant VCF index"
    }
    ref_name: {
      description: "Reference genome short name"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    sawfish_expected_bed_male: {
      description: "Expected allosome copy number BED for XY samples"
    }
    sawfish_expected_bed_female: {
      description: "Expected allosome copy number BED for XX samples"
    }
    methbat_region_tsv: {
      description: "Regions for MethBat methylation profiling in tab-separated format"
    }
    run_starphase: {
      description: "Whether to run StarPhase task"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    String sex
    File aligned_hifi_reads
    File aligned_hifi_reads_index
    File? aligned_fail_reads
    File? aligned_fail_reads_index
    File trgt_catalog
    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index
    String ref_name
    File ref_fasta
    File ref_index
    File sawfish_expected_bed_male
    File sawfish_expected_bed_female
    File methbat_region_tsv
    Boolean run_starphase
    RuntimeAttributes default_runtime_attributes
  }

  Array[File] hiphase_input_vcfs = [
    small_variant_vcf,
    sv_vcf
  ]
  Array[File] hiphase_input_vcf_indices = [
    small_variant_vcf_index,
    sv_vcf_index
  ]

  # In order to properly delocalize the VCF outputs of hiphase in cloud engines
  # we need to generate the names of the outputs and pass these to the call.
  scatter (vcf_index in range(length(hiphase_input_vcfs))) {
    String phased_vcf_name = basename(hiphase_input_vcfs[vcf_index], ".vcf.gz") + ".phased.vcf.gz"
    String phased_vcf_index_name = basename(hiphase_input_vcf_indices[vcf_index], ".vcf.gz.tbi") + ".phased.vcf.gz.tbi"
  }

  String haplotagged_bam_name = "~{sample_id}.~{ref_name}.haplotagged.bam"
  String haplotagged_bam_index_name = "~{sample_id}.~{ref_name}.haplotagged.bam.bai"

  call Hiphase.hiphase { input:
    sample_id = sample_id,
    vcfs = hiphase_input_vcfs,
    vcf_indices = hiphase_input_vcf_indices,
    phased_vcf_names = phased_vcf_name,
    phased_vcf_index_names = phased_vcf_index_name,
    aligned_bams = [
      aligned_hifi_reads
    ],
    aligned_bam_indices = [
      aligned_hifi_reads_index
    ],
    haplotagged_bam_names = [
      haplotagged_bam_name
    ],
    haplotagged_bam_index_names = [
      haplotagged_bam_index_name
    ],
    ref_name = ref_name,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    runtime_attributes = default_runtime_attributes
  }

  # hiphase.phased_vcfs[0] -> phased small variant VCF
  # hiphase.phased_vcfs[1] -> phased SV VCF
  # hiphase.haplotagged_bams[0]/haplotagged_bam_indices[0] -> the single
  # merged BAM downstream.wdl always hands hiphase

  call Trgt.trgt_genotype { input:
    sample_id = sample_id,
    sex = sex,
    aligned_bam = hiphase.haplotagged_bams[0],
    aligned_bam_index = hiphase.haplotagged_bam_indices[0],
    fail_reads_bam = aligned_fail_reads,
    fail_reads_bam_index = aligned_fail_reads_index,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    trgt_bed = trgt_catalog,
    expected_male_bed = sawfish_expected_bed_male,
    expected_female_bed = sawfish_expected_bed_female,
    out_prefix = "~{sample_id}.~{ref_name}",
    min_read_quality = -1.0,
    max_depth = 150,
    runtime_attributes = default_runtime_attributes
  }

  call Pbjam.pbjam_bam_stats { input:
    sample_id = sample_id,
    ref_name = ref_name,
    bam = hiphase.haplotagged_bams[0],
    bam_index = hiphase.haplotagged_bam_indices[0],
    runtime_attributes = default_runtime_attributes
  }

  call Bcftools.bcftools_stats_roh_small_variants { input:
    sample_id = sample_id,
    vcf = hiphase.phased_vcfs[0],
    ref_fasta = ref_fasta,
    ref_name = ref_name,
    runtime_attributes = default_runtime_attributes
  }

  call Bcftools.sv_stats { input:
    sample_id = sample_id,
    ref_name = ref_name,
    vcf = hiphase.phased_vcfs[1],
    runtime_attributes = default_runtime_attributes
  }

  call Methbat.methbat_pileup { input:
    haplotagged_bam = hiphase.haplotagged_bams[0],
    haplotagged_bam_index = hiphase.haplotagged_bam_indices[0],
    out_prefix = "~{sample_id}.~{ref_name}",
    runtime_attributes = default_runtime_attributes
  }

  if (defined(methbat_pileup.cpg_pileup_bed)) {
    # If any cpg_pileup_beds are generated, we can run methbat
    call Methbat.methbat_profile as methbat_profile_task { input:
      cpg_pileup_bed = select_first([
        methbat_pileup.cpg_pileup_bed
      ]),
      region_tsv = methbat_region_tsv,
      out_prefix = "~{sample_id}.~{ref_name}",
      runtime_attributes = default_runtime_attributes
    }
  }

  if (run_starphase) {
    call Pbstarphase.pbstarphase_diplotype { input:
      out_prefix = sample_id,
      phased_small_variant_vcf = hiphase.phased_vcfs[0],
      phased_small_variant_vcf_index = hiphase.phased_vcf_indices[0],
      phased_structural_variant_vcf = hiphase.phased_vcfs[1],
      phased_structural_variant_vcf_index = hiphase.phased_vcf_indices[1],
      aligned_bam = hiphase.haplotagged_bams[0],
      aligned_bam_index = hiphase.haplotagged_bam_indices[0],
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      runtime_attributes = default_runtime_attributes
    }
  }

  output {
    # hiphase outputs
    File merged_haplotagged_bam = hiphase.haplotagged_bams[0]
    File merged_haplotagged_bam_index = hiphase.haplotagged_bam_indices[0]
    File phased_small_variant_vcf = hiphase.phased_vcfs[0]
    File phased_small_variant_vcf_index = hiphase.phased_vcf_indices[0]
    File phased_sv_vcf = hiphase.phased_vcfs[1]
    File phased_sv_vcf_index = hiphase.phased_vcf_indices[1]
    File phase_stats = hiphase.phase_stats
    File phase_blocks = hiphase.phase_blocks
    File phase_haplotags = hiphase.phase_haplotags
    String stat_phased_basepairs = hiphase.stat_phased_basepairs
    String stat_phase_block_ng50 = hiphase.stat_phase_block_ng50

    # bam stats
    File read_length_plot = pbjam_bam_stats.read_length_plot
    File read_quality_plot = pbjam_bam_stats.read_quality_plot
    File mapq_distribution_plot = pbjam_bam_stats.mapq_distribution_plot
    File mg_distribution_plot = pbjam_bam_stats.mg_distribution_plot
    String stat_read_count = pbjam_bam_stats.stat_read_count
    String stat_read_length_mean = pbjam_bam_stats.stat_read_length_mean
    String stat_read_length_median = pbjam_bam_stats.stat_read_length_median
    String stat_read_length_n50 = pbjam_bam_stats.stat_read_length_n50
    String stat_read_quality_mean = pbjam_bam_stats.stat_read_quality_mean
    String stat_read_quality_median = pbjam_bam_stats.stat_read_quality_median
    String stat_mapped_read_count = pbjam_bam_stats.stat_mapped_read_count
    String stat_mapped_read_percent = pbjam_bam_stats.stat_mapped_read_percent
    String stat_gap_compressed_identity_mean = pbjam_bam_stats.stat_gap_compressed_identity_mean
    String stat_gap_compressed_identity_median = pbjam_bam_stats.stat_gap_compressed_identity_median
    File trgt_coverage_dropouts = trgt_genotype.dropouts

    # small variant stats
    File small_variant_stats = bcftools_stats_roh_small_variants.stats
    File bcftools_roh_out = bcftools_stats_roh_small_variants.roh_out
    File bcftools_roh_bed = bcftools_stats_roh_small_variants.roh_bed
    String stat_SNV_count = bcftools_stats_roh_small_variants.stat_SNV_count
    String stat_INDEL_count = bcftools_stats_roh_small_variants.stat_INDEL_count
    String stat_TSTV_ratio = bcftools_stats_roh_small_variants.stat_TSTV_ratio
    String stat_HETHOM_ratio = bcftools_stats_roh_small_variants.stat_HETHOM_ratio
    File snv_distribution_plot = bcftools_stats_roh_small_variants.snv_distribution_plot
    File indel_distribution_plot = bcftools_stats_roh_small_variants.indel_distribution_plot

    # sv stats
    String stat_sv_DUP_count = sv_stats.stat_sv_DUP_count
    String stat_sv_DEL_count = sv_stats.stat_sv_DEL_count
    String stat_sv_INS_count = sv_stats.stat_sv_INS_count
    String stat_sv_INV_count = sv_stats.stat_sv_INV_count
    String stat_sv_BND_count = sv_stats.stat_sv_BND_count
    String stat_sv_SWAP_count = sv_stats.stat_sv_SWAP_count
    File sv_stats_plot = sv_stats.sv_stats_plot

    # trgt outputs
    File trgt_vcf = trgt_genotype.vcf
    File trgt_vcf_index = trgt_genotype.vcf_index
    File trgt_spanning_reads = trgt_genotype.bam
    File trgt_spanning_reads_index = trgt_genotype.bam_index
    String stat_trgt_genotyped_count = trgt_genotype.stat_genotyped_count
    String stat_trgt_uncalled_count = trgt_genotype.stat_uncalled_count

    # methylation outputs and profile
    File? cpg_pileup_bed = methbat_pileup.cpg_pileup_bed
    File? cpg_pileup_bed_index = methbat_pileup.cpg_pileup_bed_index
    File? hmcpg_pileup_bed = methbat_pileup.hmcpg_pileup_bed
    File? hmcpg_pileup_bed_index = methbat_pileup.hmcpg_pileup_bed_index
    String stat_cpg_hap1_count = methbat_pileup.stat_hap1_cpg_count
    String stat_cpg_hap2_count = methbat_pileup.stat_hap2_cpg_count
    String stat_cpg_combined_count = methbat_pileup.stat_combined_cpg_count
    File? methbat_profile = methbat_profile_task.profile
    String stat_methbat_methylated_count = select_first([
      methbat_profile_task.stat_methbat_methylated_count,
      "0"
    ])
    String stat_methbat_unmethylated_count = select_first([
      methbat_profile_task.stat_methbat_unmethylated_count,
      "0"
    ])
    String stat_methbat_asm_count = select_first([
      methbat_profile_task.stat_methbat_asm_count,
      "0"
    ])

    # pbstarphase outputs
    File? pbstarphase_json = pbstarphase_diplotype.summary_json
    File? pbstarphase_tsv = pbstarphase_diplotype.pharmcat_tsv

    # qc messages
    Array[String] msg = flatten([
      trgt_genotype.msg
    ])
  }
}

