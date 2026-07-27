version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "../wdl-common/wdl/tasks/glnexus.wdl" as Glnexus
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish

workflow joint {
  meta {
    description: "Tasks for joint-calling variants from a set of samples and splitting the joint calls by sample for parallel phasing."
    outputs: {
      split_joint_structural_variant_vcfs: {
        description: "Joint-called structural variant VCFs, split by sample"
      },
      split_joint_structural_variant_vcf_indices: {
        description: "Joint-called structural variant VCF indices, split by sample"
      },
      split_joint_small_variant_vcfs: {
        description: "Joint-called small variant VCFs, split by sample"
      },
      split_joint_small_variant_vcf_indices: {
        description: "Joint-called small variant VCF indices, split by sample"
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
      }
    }
  }

  parameter_meta {
    family_id: {
      description: "Cohort ID"
    }
    sample_ids: {
      description: "Sample IDs"
    }
    gvcfs: {
      description: "GVCFs"
    }
    gvcf_indices: {
      description: "GVCF indices"
    }
    discover_tars: {
      description: "Sawfish discover output tarballs"
    }
    aligned_bams: {
      description: "Aligned BAMs"
    }
    aligned_bam_indices: {
      description: "Aligned BAM indices"
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
    glnexus_mem_gb: {
      description: "GLnexus memory (GB)"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String family_id
    Array[String] sample_ids
    Array[File] gvcfs
    Array[File] gvcf_indices
    Array[File] discover_tars
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    String ref_name
    File ref_fasta
    File ref_index
    Int glnexus_mem_gb
    RuntimeAttributes default_runtime_attributes
  }

  # In order to properly delocalize the outputs of sawfish_call in cloud engines
  # we need to generate the names of the outputs and pass these to the call.
  scatter (sample_id in sample_ids) {
    String copynum_bedgraph_name = "~{sample_id}.~{family_id}.joint.~{ref_name}.structural_variants.copynum.bedgraph"
    String depth_bw_name = "~{sample_id}.~{family_id}.joint.~{ref_name}.structural_variants.depth.bw"
    String gc_bias_corrected_depth_bw_name = "~{sample_id}.~{family_id}.joint.~{ref_name}.structural_variants.gc_bias_corrected_depth.bw"
    String copynum_summary_name = "~{sample_id}.~{family_id}.joint.~{ref_name}.structural_variants.copynum.summary.json"
  }

  call Sawfish.sawfish_call { input:
    sample_ids = sample_ids,
    discover_tars = discover_tars,
    aligned_bams = aligned_bams,
    aligned_bam_indices = aligned_bam_indices,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    out_prefix = "~{family_id}.joint.~{ref_name}.structural_variants",
    copynum_bedgraph_names = copynum_bedgraph_name,
    depth_bw_names = depth_bw_name,
    gc_bias_corrected_depth_bw_names = gc_bias_corrected_depth_bw_name,
    copynum_summary_names = copynum_summary_name,
    runtime_attributes = default_runtime_attributes
  }

  # In order to properly delocalize the outputs of split_vcf_by_sample in cloud engines
  # we need to generate the names of the outputs and pass these to the call.
  String sv_vcf_basename = basename(sawfish_call.vcf, ".vcf.gz")
  scatter (sample_id in sample_ids) {
    String split_sv_vcf_name = "~{sample_id}.~{sv_vcf_basename}.vcf.gz"
    String split_sv_vcf_index_name = "~{sample_id}.~{sv_vcf_basename}.vcf.gz.tbi"
  }

  call Bcftools.split_vcf_by_sample as split_sawfish { input:
    sample_ids = sample_ids,
    vcf = sawfish_call.vcf,
    vcf_index = sawfish_call.vcf_index,
    split_vcf_names = split_sv_vcf_name,
    split_vcf_index_names = split_sv_vcf_index_name,
    exclude_uncalled = false,
    runtime_attributes = default_runtime_attributes
  }

  call Glnexus.glnexus { input:
    cohort_id = family_id + ".joint",
    gvcfs = gvcfs,
    gvcf_indices = gvcf_indices,
    ref_name = ref_name,
    mem_gb = glnexus_mem_gb,
    runtime_attributes = default_runtime_attributes
  }

  # In order to properly delocalize the outputs of split_vcf_by_sample in cloud engines
  # we need to generate the names of the outputs and pass these to the call.
  String glnexus_vcf_basename = basename(glnexus.vcf, ".vcf.gz")
  scatter (sample_id in sample_ids) {
    String split_glnexus_vcf_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz"
    String split_glnexus_vcf_index_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz.tbi"
  }

  call Bcftools.split_vcf_by_sample as split_glnexus { input:
    sample_ids = sample_ids,
    vcf = glnexus.vcf,
    vcf_index = glnexus.vcf_index,
    split_vcf_names = split_glnexus_vcf_name,
    split_vcf_index_names = split_glnexus_vcf_index_name,
    runtime_attributes = default_runtime_attributes
  }

  output {
    Array[File] split_joint_structural_variant_vcfs = split_sawfish.split_vcfs
    Array[File] split_joint_structural_variant_vcf_indices = split_sawfish.split_vcf_indices
    Array[File] split_joint_small_variant_vcfs = split_glnexus.split_vcfs
    Array[File] split_joint_small_variant_vcf_indices = split_glnexus.split_vcf_indices
    File sv_supporting_reads = select_first([
      sawfish_call.supporting_reads
    ])
    Array[File] sv_copynum_bedgraph = sawfish_call.copynum_bedgraph
    Array[File] sv_depth_bw = sawfish_call.depth_bw
    Array[File] sv_gc_bias_corrected_depth_bw = sawfish_call.gc_bias_corrected_depth_bw
    Array[File] sv_copynum_summary = sawfish_call.copynum_summary
  }
}

