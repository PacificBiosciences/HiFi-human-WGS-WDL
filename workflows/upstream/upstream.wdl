version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/workflows/pbmm2/pbmm2.wdl" as Pbmm2
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../wdl-common/wdl/workflows/deepvariant/deepvariant.wdl" as DeepVariant
import "../wdl-common/wdl/tasks/samtools.wdl" as Samtools
import "../wdl-common/wdl/tasks/mosdepth.wdl" as Mosdepth
import "../wdl-common/wdl/tasks/paraphase.wdl" as Paraphase
import "../wdl-common/wdl/tasks/mitorsaw.wdl" as Mitorsaw

workflow upstream {
  meta {
    description: "Given a set of HiFi reads for a human sample, run steps upstream of phasing."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    hifi_reads: {
      name: "HiFi reads (BAMs)"
    }
    fail_reads: {
      name: "Failed reads (BAMs)"
    }
    fail_reads_bed: {
      name: "Subset of genome for which to include fail reads; BED"
    }
    fail_reads_bait_fasta: {
      name: "FASTA of reference sequences for baiting fail reads; FASTA"
    }
    fail_reads_bait_index: {
      name: "Index of reference sequences for baiting fail reads; FASTA index"
    }
    ref_map_file: {
      name: "TSV containing reference genome information"
    }
    max_reads_per_alignment_chunk: {
      name: "Maximum reads per alignment chunk"
    }
    single_sample: {
      name: "Single sample workflow"
    }
    gpu: {
      name: "Use GPU for DeepVariant"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String? sex
    Array[File] hifi_reads
    Array[File]? fail_reads

    File? fail_reads_bed
    File? fail_reads_bait_fasta
    File? fail_reads_bait_index

    File ref_map_file

    Int max_reads_per_alignment_chunk

    Boolean single_sample = false

    Boolean gpu

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)

  if (
    defined(fail_reads) && defined(fail_reads_bed) && defined(fail_reads_bait_fasta) && defined(fail_reads_bait_index)
  ) {
    scatter (fail_read_bam in select_first([fail_reads])) {
      call Pbmm2.pbmm2_align_wgs as bait_fail_reads {
        input:
          sample_id          = sample_id,
          bam                = fail_read_bam,
          ref_fasta          = select_first([fail_reads_bait_fasta]),
          ref_index          = select_first([fail_reads_bait_index]),
          ref_name           = 'bait_ref',
          keep_unmapped      = false,
          min_length         = 1000,
          runtime_attributes = default_runtime_attributes
      }

      call Pbmm2.pbmm2_align_wgs as align_captured_fail_reads {
        input:
          sample_id          = sample_id,
          bam                = bait_fail_reads.aligned_bam,
          ref_fasta          = ref_map["fasta"],            # !FileCoercion
          ref_index          = ref_map["fasta_index"],      # !FileCoercion
          ref_name           = ref_map["name"],
          keep_unmapped      = false,
          runtime_attributes = default_runtime_attributes
      }

      call Samtools.subset_bam {
        input:
          bed                = select_first([fail_reads_bed]),
          aligned_bam        = align_captured_fail_reads.aligned_bam,
          aligned_bam_index  = align_captured_fail_reads.aligned_bam_index,
          ref_index          = ref_map["fasta_index"],                                 # !FileCoercion
          out_prefix         = "~{sample_id}.~{ref_map['name']}.baited_fail_reads",
          runtime_attributes = default_runtime_attributes
      }
    }

    # merge aligned bams if there are multiple
    if (length(subset_bam.bam) > 1) {
      call Samtools.samtools_merge as merge_fail_reads {
        input:
          bams               = subset_bam.bam,
          out_prefix         = "~{sample_id}.~{ref_map['name']}.fail_reads",
          runtime_attributes = default_runtime_attributes
      }
    }

    # select the merged bam if it exists, otherwise select the first (only) aligned bam
    File aligned_fail_bam_data  = select_first([merge_fail_reads.merged_bam, subset_bam.bam[0]])
    File aligned_fail_bam_index = select_first([merge_fail_reads.merged_bam_index, subset_bam.bam[0]])
  }

  String include_fail_reads = 
    if (defined(aligned_fail_bam_data)) 
    then "Including fail_reads for TRGT genotyping for regions specified in the TRGT catalog."
    else ""

  scatter (hifi_read_bam in hifi_reads) {
    call Pbmm2.pbmm2 as pbmm2 {
      input:
        sample_id                  = sample_id,
        bam                        = hifi_read_bam,
        max_reads_per_chunk        = max_reads_per_alignment_chunk,
        ref_fasta                  = ref_map["fasta"],              # !FileCoercion
        ref_index                  = ref_map["fasta_index"],        # !FileCoercion
        ref_name                   = ref_map["name"],
        default_runtime_attributes = default_runtime_attributes
    }
  }

  # merge aligned bams if there are multiple
  if (length(flatten(pbmm2.aligned_bams)) > 1) {
    call Samtools.samtools_merge as merge_hifi_reads {
      input:
        bams               = flatten(pbmm2.aligned_bams),
        out_prefix         = "~{sample_id}.~{ref_map['name']}.hifi_reads",
        runtime_attributes = default_runtime_attributes
    }
  }

  # select the merged bam if it exists, otherwise select the first (only) aligned bam
  File aligned_bam_data  = select_first([merge_hifi_reads.merged_bam, flatten(pbmm2.aligned_bams)[0]])
  File aligned_bam_index = select_first([merge_hifi_reads.merged_bam_index, flatten(pbmm2.aligned_bam_indices)[0]])

  call Mosdepth.mosdepth {
    input:
      sample_id          = sample_id,
      ref_name           = ref_map["name"],
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      infer_sex          = true,
      runtime_attributes = default_runtime_attributes
  }

  String qc_sex = 
    if (defined(sex) && (mosdepth.inferred_sex != sex)) 
    then "~{sample_id}: Reported sex ~{sex} does not match inferred sex ~{mosdepth.inferred_sex}."
    else ""

  call DeepVariant.deepvariant {
    input:
      sample_id                  = sample_id,
      aligned_bams               = [aligned_bam_data],
      aligned_bam_indices        = [aligned_bam_index],
      ref_fasta                  = ref_map["fasta"],       # !FileCoercion
      ref_index                  = ref_map["fasta_index"], # !FileCoercion
      ref_name                   = ref_map["name"],
      gpu                        = gpu,
      default_runtime_attributes = default_runtime_attributes
  }

  call Sawfish.sawfish_discover {
    input:
      sample_id               = sample_id,
      sex                     = mosdepth.inferred_sex,
      aligned_bam             = aligned_bam_data,
      aligned_bam_index       = aligned_bam_index,
      ref_fasta               = ref_map["fasta"],                       # !FileCoercion
      ref_index               = ref_map["fasta_index"],                 # !FileCoercion
      exclude_bed             = ref_map["sawfish_exclude_bed"],         # !FileCoercion
      exclude_bed_index       = ref_map["sawfish_exclude_bed_index"],   # !FileCoercion
      expected_male_bed       = ref_map["sawfish_expected_bed_male"],   # !FileCoercion
      expected_female_bed     = ref_map["sawfish_expected_bed_female"], # !FileCoercion
      small_variant_vcf       = deepvariant.vcf,
      small_variant_vcf_index = deepvariant.vcf_index,
      out_prefix              = "~{sample_id}",
      runtime_attributes      = default_runtime_attributes
  }

  call Paraphase.paraphase {
    input:
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],         # !FileCoercion
      ref_index          = ref_map["fasta_index"],   # !FileCoercion
      sample_id          = sample_id,
      runtime_attributes = default_runtime_attributes
  }

  call Mitorsaw.mitorsaw {
    input:
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],                  # !FileCoercion
      ref_index          = ref_map["fasta_index"],            # !FileCoercion
      out_prefix         = "~{sample_id}.~{ref_map['name']}",
      runtime_attributes = default_runtime_attributes
  }

  if (single_sample) {
    String copynum_bedgraph_name           = "~{sample_id}.~{ref_map['name']}.structural_variants.copynum.bedgraph"
    String depth_bw_name                   = "~{sample_id}.~{ref_map['name']}.structural_variants.depth.bw"
    String gc_bias_corrected_depth_bw_name = "~{sample_id}.~{ref_map['name']}.structural_variants.gc_bias_corrected_depth.bw"
    String maf_bw_name                     = "~{sample_id}.~{ref_map['name']}.structural_variants.maf.bw"
    String copynum_summary_name            = "~{sample_id}.~{ref_map['name']}.structural_variants.copynum.summary.json"


    call Sawfish.sawfish_call {
      input:
        sample_ids                       = [sample_id],
        discover_tars                    = [sawfish_discover.discover_tar],
        aligned_bams                     = [aligned_bam_data],
        aligned_bam_indices              = [aligned_bam_index],
        ref_fasta                        = ref_map["fasta"],                                      # !FileCoercion
        ref_index                        = ref_map["fasta_index"],                                # !FileCoercion
        out_prefix                       = "~{sample_id}.~{ref_map['name']}.structural_variants",
        copynum_bedgraph_names           = [copynum_bedgraph_name],
        depth_bw_names                   = [depth_bw_name],
        gc_bias_corrected_depth_bw_names = [gc_bias_corrected_depth_bw_name],
        maf_bw_names                     = [maf_bw_name],
        copynum_summary_names            = [copynum_summary_name],
        runtime_attributes               = default_runtime_attributes
    }

    File copynum_bedgraph_output           = sawfish_call.copynum_bedgraph[0]
    File depth_bw_output                   = sawfish_call.depth_bw[0]
    File gc_bias_corrected_depth_bw_output = sawfish_call.gc_bias_corrected_depth_bw[0]
    File maf_bw_output                     = sawfish_call.maf_bw[0]
    File copynum_summary_output            = sawfish_call.copynum_summary[0]
  }

  output {
    # alignments
    File  aligned_hifi_reads       = aligned_bam_data
    File  aligned_hifi_reads_index = aligned_bam_index
    File? aligned_fail_reads       = aligned_fail_bam_data
    File? aligned_fail_reads_index = aligned_fail_bam_index

    # mosdepth outputs
    File   mosdepth_summary                 = mosdepth.summary
    File   mosdepth_region_bed              = mosdepth.region_bed
    File   mosdepth_region_bed_index        = mosdepth.region_bed_index
    File   mosdepth_depth_distribution_plot = mosdepth.depth_distribution_plot
    String inferred_sex                     = mosdepth.inferred_sex
    String stat_depth_mean                  = mosdepth.stat_depth_mean

    # per sample sv signatures
    File discover_tar = sawfish_discover.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf                        = sawfish_call.vcf
    File? sv_vcf_index                  = sawfish_call.vcf_index
    File? sv_supporting_reads           = sawfish_call.supporting_reads
    File? sv_copynum_bedgraph           = copynum_bedgraph_output
    File? sv_depth_bw                   = depth_bw_output
    File? sv_gc_bias_corrected_depth_bw = gc_bias_corrected_depth_bw_output
    File? sv_maf_bw                     = maf_bw_output
    File? sv_copynum_summary            = copynum_summary_output

    # small variant outputs
    File small_variant_vcf        = deepvariant.vcf
    File small_variant_vcf_index  = deepvariant.vcf_index
    File small_variant_gvcf       = deepvariant.gvcf
    File small_variant_gvcf_index = deepvariant.gvcf_index

    # paraphase outputs
    File? paraphase_output_json         = paraphase.out_json
    File? paraphase_realigned_bam       = paraphase.bam
    File? paraphase_realigned_bam_index = paraphase.bam_index
    File? paraphase_vcfs                = paraphase.vcfs_tar

    # per sample mitorsaw outputs
    File mitorsaw_vcf       = mitorsaw.vcf
    File mitorsaw_vcf_index = mitorsaw.vcf_index
    File mitorsaw_hap_stats = mitorsaw.hap_stats

    # qc messages
    Array[String] msg = flatten(
      [
        flatten(pbmm2.msg),
        [include_fail_reads],
        [qc_sex],
        sawfish_discover.msg
      ]
    )
  }
}
