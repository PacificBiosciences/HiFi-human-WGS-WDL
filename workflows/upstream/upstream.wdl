version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/kivvi.wdl" as Kivvi
import "../wdl-common/wdl/tasks/mitorsaw.wdl" as Mitorsaw
import "../wdl-common/wdl/tasks/mosdepth.wdl" as Mosdepth
import "../wdl-common/wdl/tasks/paraphase.wdl" as Paraphase
import "../wdl-common/wdl/tasks/pbsamoa.wdl" as Pbsamoa
import "../wdl-common/wdl/tasks/samtools.wdl" as Samtools
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../wdl-common/wdl/workflows/deepvariant/deepvariant.wdl" as DeepVariant
import "../wdl-common/wdl/workflows/parabricks/parabricks_deepvariant.wdl" as ParabricksDeepVariant
import "../wdl-common/wdl/workflows/pbmm2/pbmm2.wdl" as Pbmm2

workflow upstream {
  meta {
    description: "Align reads and call variants."
    help: "Called by singleton and family entrypoints.  Not intended to be called directly."
    category: "Subworkflow"
    outputs: {
      aligned_hifi_reads: {
        description: "Aligned HiFi reads"
      },
      aligned_hifi_reads_index: {
        description: "Index for aligned HiFi reads"
      },
      aligned_fail_reads: {
        description: "Aligned fail reads"
      },
      aligned_fail_reads_index: {
        description: "Index for aligned fail reads"
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
      inferred_sex: {
        description: "Inferred sex"
      },
      stat_depth_mean: {
        description: "Mean depth"
      },
      discover_tar: {
        description: "Sawfish discover outputs"
      },
      sv_vcf: {
        description: "Structural variant VCF"
      },
      sv_vcf_index: {
        description: "Index for structural variant VCF"
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
      small_variant_vcf: {
        description: "Small variant VCF"
      },
      small_variant_vcf_index: {
        description: "Index for small variant VCF"
      },
      small_variant_gvcf: {
        description: "Small variant GVCF"
      },
      small_variant_gvcf_index: {
        description: "Index for small variant GVCF"
      },
      paraphase_output_json: {
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
      msg: {
        description: "Messages from the workflow"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    hifi_reads: {
      description: "Unaligned hifi_reads BAMs"
    }
    fail_reads: {
      description: "Unaligned fail_reads BAMs",
      group: "Common"
    }
    fail_reads_bed: {
      description: "Subset of genome for which to include fail reads",
      group: "Common"
    }
    fail_reads_bait_index: {
      description: "Index of reference sequences for baiting fail reads",
      group: "Common"
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
    max_norm_female_chrY_depth: {
      description: "Maximum expected normalized chrY depth for samples without chrY"
    }
    paraphase_genome_build: {
      description: "Genome reference build parameter for Paraphase"
    }
    sawfish_exclude_bed: {
      description: "Regions to be excluded for Sawfish CNV calls in gzipped BED format"
    }
    sawfish_exclude_bed_index: {
      description: "Regions to be excluded for Sawfish CNV calls in gzipped BED format index file"
    }
    sawfish_expected_bed_male: {
      description: "Expected allosome copy number BED for XY samples"
    }
    sawfish_expected_bed_female: {
      description: "Expected allosome copy number BED for XX samples"
    }
    use_alignment_chunking: {
      description: "Whether to chunk BAM files for alignment. If false, all reads will be aligned in a single chunk."
    }
    single_sample: {
      description: "Single sample workflow",
      group: "Common"
    }
    use_gpu: {
      description: "Use GPU when possible"
    }
    use_parabricks_deepvariant: {
      description: "Use Parabricks DeepVariant"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] hifi_reads
    Array[File]? fail_reads
    File? fail_reads_bed
    File? fail_reads_bait_index
    String ref_name
    File ref_fasta
    File ref_index
    Float max_norm_female_chrY_depth
    String paraphase_genome_build
    File sawfish_exclude_bed
    File sawfish_exclude_bed_index
    File sawfish_expected_bed_male
    File sawfish_expected_bed_female
    Boolean use_alignment_chunking
    Boolean single_sample = false
    Boolean use_gpu
    Boolean use_parabricks_deepvariant
    RuntimeAttributes default_runtime_attributes
  }

  scatter (hifi_read_bam in hifi_reads) {
    call Pbmm2.pbmm2 as pbmm2 { input:
      sample_id = sample_id,
      bam = hifi_read_bam,
      use_alignment_chunking = use_alignment_chunking,
      ref_fasta = ref_fasta,
      ref_name = ref_name,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  # merge aligned bams if there are multiple
  if (length(flatten(pbmm2.aligned_bams)) > 1) {
    call Pbsamoa.pbsamoa_merge as merge_hifi_reads { input:
      bams = flatten(pbmm2.aligned_bams),
      out_prefix = "~{sample_id}.~{ref_name}.hifi_reads",
      runtime_attributes = default_runtime_attributes
    }
  }

  # select the merged bam if it exists, otherwise select the first (only) aligned bam
  File aligned_bam_data = select_first([
    merge_hifi_reads.merged_bam,
    flatten(pbmm2.aligned_bams)[0]
  ])
  File aligned_bam_index = select_first([
    merge_hifi_reads.merged_bam_index,
    flatten(pbmm2.aligned_bam_indices)[0]
  ])

  call Mosdepth.mosdepth { input:
    sample_id = sample_id,
    ref_name = ref_name,
    aligned_bam = aligned_bam_data,
    aligned_bam_index = aligned_bam_index,
    infer_sex = true,
    max_norm_female_chrY_depth = max_norm_female_chrY_depth,
    runtime_attributes = default_runtime_attributes
  }

  if (use_gpu && use_parabricks_deepvariant) {
    call ParabricksDeepVariant.parabricks_deepvariant { input:
      sample_id = sample_id,
      aligned_bam = aligned_bam_data,
      aligned_bam_index = aligned_bam_index,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_name = ref_name,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  if (!use_parabricks_deepvariant || (!use_gpu && use_parabricks_deepvariant)) {
    call DeepVariant.deepvariant { input:
      sample_id = sample_id,
      aligned_bams = [
        aligned_bam_data
      ],
      aligned_bam_indices = [
        aligned_bam_index
      ],
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      ref_name = ref_name,
      gpu = use_gpu,
      default_runtime_attributes = default_runtime_attributes
    }
  }

  File deepvariant_vcf = select_first([
    parabricks_deepvariant.vcf,
    deepvariant.vcf
  ])
  File deepvariant_vcf_index = select_first([
    parabricks_deepvariant.vcf_index,
    deepvariant.vcf_index
  ])
  File deepvariant_gvcf = select_first([
    parabricks_deepvariant.gvcf,
    deepvariant.gvcf
  ])
  File deepvariant_gvcf_index = select_first([
    parabricks_deepvariant.gvcf_index,
    deepvariant.gvcf_index
  ])

  call Sawfish.sawfish_discover { input:
    sample_id = sample_id,
    sex = mosdepth.inferred_sex,
    aligned_bam = aligned_bam_data,
    aligned_bam_index = aligned_bam_index,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    exclude_bed = sawfish_exclude_bed,
    exclude_bed_index = sawfish_exclude_bed_index,
    expected_male_bed = sawfish_expected_bed_male,
    expected_female_bed = sawfish_expected_bed_female,
    out_prefix = sample_id,
    runtime_attributes = default_runtime_attributes
  }

  if (defined(fail_reads) && length(select_first([
    fail_reads
  ])) > 0 && defined(fail_reads_bed) && defined(fail_reads_bait_index)) {
    scatter (fail_read_bam in select_first([
      fail_reads
    ])) {
      call Pbmm2.pbmm2_align_wgs as bait_fail_reads { input:
        sample_id = sample_id,
        bam = fail_read_bam,
        pbmm2_index = select_first([
          fail_reads_bait_index
        ]),
        ref_name = "bait_ref",
        keep_unmapped = false,
        min_length = 1000,
        runtime_attributes = default_runtime_attributes
      }

      call Pbmm2.create_pbmm2_index { input:
        ref_fasta = ref_fasta,
        runtime_attributes = default_runtime_attributes
      }

      call Pbmm2.pbmm2_align_wgs as align_captured_fail_reads { input:
        sample_id = sample_id,
        bam = bait_fail_reads.aligned_bam,
        pbmm2_index = create_pbmm2_index.index,
        ref_name = ref_name,
        keep_unmapped = false,
        runtime_attributes = default_runtime_attributes
      }

      call Samtools.subset_bam { input:
        bed = select_first([
          fail_reads_bed
        ]),
        aligned_bam = align_captured_fail_reads.aligned_bam,
        aligned_bam_index = align_captured_fail_reads.aligned_bam_index,
        ref_index = ref_index,
        out_prefix = "~{sample_id}.~{ref_name}.baited_fail_reads",
        runtime_attributes = default_runtime_attributes
      }
    }

    # merge aligned bams if there are multiple
    if (length(subset_bam.bam) > 1) {
      call Pbsamoa.pbsamoa_merge as merge_fail_reads { input:
        bams = subset_bam.bam,
        out_prefix = "~{sample_id}.~{ref_name}.fail_reads",
        runtime_attributes = default_runtime_attributes
      }
    }

    # select the merged bam if it exists, otherwise select the first (only) aligned bam
    File aligned_fail_bam_data = select_first([
      merge_fail_reads.merged_bam,
      subset_bam.bam[0]
    ])
    File aligned_fail_bam_index = select_first([
      merge_fail_reads.merged_bam_index,
      subset_bam.bam_index[0]
    ])
  }

  String include_fail_reads = if (defined(aligned_fail_bam_data))
    then "Including fail_reads for TRGT genotyping for regions specified in the TRGT catalog."
    else ""

  call Paraphase.paraphase { input:
    aligned_bam = aligned_bam_data,
    aligned_bam_index = aligned_bam_index,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    genome = paraphase_genome_build,
    sample_id = sample_id,
    runtime_attributes = default_runtime_attributes
  }

  call Mitorsaw.mitorsaw { input:
    aligned_bam = aligned_bam_data,
    aligned_bam_index = aligned_bam_index,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    out_prefix = "~{sample_id}.~{ref_name}",
    runtime_attributes = default_runtime_attributes
  }

  call Kivvi.kivvi_kiv2 { input:
    aligned_bam = aligned_bam_data,
    aligned_bam_index = aligned_bam_index,
    out_prefix = "~{sample_id}.~{ref_name}",
    runtime_attributes = default_runtime_attributes
  }

  call Kivvi.kivvi_d4z4 { input:
    aligned_bam = aligned_bam_data,
    aligned_bam_index = aligned_bam_index,
    out_prefix = "~{sample_id}.~{ref_name}",
    runtime_attributes = default_runtime_attributes
  }

  if (single_sample) {
    String copynum_bedgraph_name = "~{sample_id}.~{ref_name}.structural_variants.copynum.bedgraph"
    String depth_bw_name = "~{sample_id}.~{ref_name}.structural_variants.depth.bw"
    String gc_bias_corrected_depth_bw_name = "~{sample_id}.~{ref_name}.structural_variants.gc_bias_corrected_depth.bw"
    String copynum_summary_name = "~{sample_id}.~{ref_name}.structural_variants.copynum.summary.json"

    call Sawfish.sawfish_call { input:
      sample_ids = [
        sample_id
      ],
      discover_tars = [
        sawfish_discover.discover_tar
      ],
      aligned_bams = [
        aligned_bam_data
      ],
      aligned_bam_indices = [
        aligned_bam_index
      ],
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      out_prefix = "~{sample_id}.~{ref_name}.structural_variants",
      copynum_bedgraph_names = [
        copynum_bedgraph_name
      ],
      depth_bw_names = [
        depth_bw_name
      ],
      gc_bias_corrected_depth_bw_names = [
        gc_bias_corrected_depth_bw_name
      ],
      copynum_summary_names = [
        copynum_summary_name
      ],
      runtime_attributes = default_runtime_attributes
    }

    File copynum_bedgraph_output = sawfish_call.copynum_bedgraph[0]
    File depth_bw_output = sawfish_call.depth_bw[0]
    File gc_bias_corrected_depth_bw_output = sawfish_call.gc_bias_corrected_depth_bw[0]
    File copynum_summary_output = sawfish_call.copynum_summary[0]
  }

  output {
    # alignments
    File aligned_hifi_reads = aligned_bam_data
    File aligned_hifi_reads_index = aligned_bam_index
    File? aligned_fail_reads = aligned_fail_bam_data
    File? aligned_fail_reads_index = aligned_fail_bam_index

    # mosdepth outputs
    File mosdepth_summary = mosdepth.summary
    File mosdepth_region_bed = mosdepth.region_bed
    File mosdepth_region_bed_index = mosdepth.region_bed_index
    File mosdepth_depth_distribution_plot = mosdepth.depth_distribution_plot
    String inferred_sex = mosdepth.inferred_sex
    String stat_depth_mean = mosdepth.stat_depth_mean

    # per sample sv signatures
    File discover_tar = sawfish_discover.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf = sawfish_call.vcf
    File? sv_vcf_index = sawfish_call.vcf_index
    File? sv_supporting_reads = sawfish_call.supporting_reads
    File? sv_copynum_bedgraph = copynum_bedgraph_output
    File? sv_depth_bw = depth_bw_output
    File? sv_gc_bias_corrected_depth_bw = gc_bias_corrected_depth_bw_output
    File? sv_copynum_summary = copynum_summary_output

    # small variant outputs
    File small_variant_vcf = deepvariant_vcf
    File small_variant_vcf_index = deepvariant_vcf_index
    File? small_variant_gvcf = deepvariant_gvcf  # !UnnecessaryQuantifier
    File? small_variant_gvcf_index = deepvariant_gvcf_index  # !UnnecessaryQuantifier

    # paraphase outputs
    File? paraphase_output_json = paraphase.summary_json
    File? paraphase_realigned_bam = paraphase.bam
    File? paraphase_realigned_bam_index = paraphase.bam_index
    File? paraphase_vcfs = paraphase.vcfs_tar

    # per sample mitorsaw outputs
    File mitorsaw_vcf = mitorsaw.vcf
    File mitorsaw_vcf_index = mitorsaw.vcf_index
    File mitorsaw_hap_stats = mitorsaw.hap_stats

    # per sample kivvi kiv2 outputs
    File? kivvi_kiv2_vcf = kivvi_kiv2.vcf
    File? kivvi_kiv2_vcf_index = kivvi_kiv2.vcf_index
    File? kivvi_kiv2_json = kivvi_kiv2.json
    File? kivvi_kiv2_realigned_bam = kivvi_kiv2.realigned_bam
    File? kivvi_kiv2_realigned_bam_index = kivvi_kiv2.realigned_bam_index
    File? kivvi_kiv2_allele_plot = kivvi_kiv2.allele_plot

    # per sample kivvi d4z4 outputs
    File? kivvi_d4z4_vcf = kivvi_d4z4.vcf
    File? kivvi_d4z4_vcf_index = kivvi_d4z4.vcf_index
    File? kivvi_d4z4_json = kivvi_d4z4.json
    File? kivvi_d4z4_realigned_bam = kivvi_d4z4.realigned_bam
    File? kivvi_d4z4_realigned_bam_index = kivvi_d4z4.realigned_bam_index
    File? kivvi_d4z4_allele_plot = kivvi_d4z4.allele_plot

    # qc messages
    Array[String] msg = flatten([
      flatten(pbmm2.msg),
      [
        include_fail_reads
      ],
      kivvi_kiv2.msg,
      kivvi_d4z4.msg,
      sawfish_discover.msg
    ])
  }
}

