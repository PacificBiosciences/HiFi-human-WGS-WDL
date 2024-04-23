version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/glnexus.wdl" as Glnexus
import "../wdl-common/wdl/tasks/pbsv.wdl" as Pbsv
import "../wdl-common/wdl/tasks/bcftools.wdl" as Bcftools


workflow joint {
  meta {
    description: "Tasks for joint-calling variants from a set of samples and splitting the joint calls by sample for parallel phasing."
  }

  parameter_meta {
    family_id: {
      name: "Cohort ID"
    }
    sample_ids: {
      name: "Sample IDs"
    }
    gvcfs: {
      name: "GVCFs"
    }   
    gvcf_indices: {
      name: "GVCF Indices"
    }
    svsigs: {
      name: "SV Signatures"
    }
    ref_map_file: {
      name: "Reference Map File"
    }
    glnexus_mem_gb: {
      name: "GLnexus Memory (GB)"
    }
    pbsv_call_mem_gb: {
      name: "PBSV Call Memory (GB)"
    }
    default_runtime_attributes: {
      name: "Default Runtime Attribute Struct"
    }
    joint_split_pbsv_vcfs: {
      name: "Joint-call PBSV VCF, split by sample"
    }
    joint_split_pbsv_vcf_indices: {
      name: "Joint-call PBSV VCF indices, split by sample"
    }
    joint_split_glenxus_vcfs: {
      name: "Joint-call GLnexus VCF, split by sample"
    }
    joint_split_glenxus_vcf_indices: {
      name: "Joint-call GLnexus VCF indices, split by sample"
    }
  }

  input {
    String family_id
    Array[String] sample_ids

    Array[File] gvcfs
    Array[File] gvcf_indices

    Array[File] svsigs

    File ref_map_file

    Int? glnexus_mem_gb
    Int? pbsv_call_mem_gb

    RuntimeAttributes default_runtime_attributes
  }

  Map  [String, String] ref_map    = read_map(ref_map_file)
  Array[Array[String]] pbsv_splits = read_json(ref_map["pbsv_splits"])  # !FileCoercion

  scatter (shard_index in range(length(pbsv_splits))) {
    Array[String] region_set = pbsv_splits[shard_index]

    call Pbsv.pbsv_call {
      input:
        sample_id          = family_id + ".joint",
        svsigs             = svsigs,
        sample_count       = length(sample_ids),
        ref_fasta          = ref_map["fasta"],       # !FileCoercion
        ref_index          = ref_map["fasta_index"], # !FileCoercion
        ref_name           = ref_map["name"],
        shard_index        = shard_index,
        regions            = region_set,
        mem_gb             = pbsv_call_mem_gb,
        runtime_attributes = default_runtime_attributes
    }
  }

    # concatenate pbsv vcfs
  call Bcftools.concat_pbsv_vcf {
    input:
      vcfs               = pbsv_call.vcf,
      vcf_indices        = pbsv_call.vcf_index,
      out_prefix         = "~{family_id}.joint.~{ref_map['name']}.structural_variants",
      runtime_attributes = default_runtime_attributes
  }

  call Bcftools.split_vcf_by_sample as split_pbsv {
    input:
      sample_ids = sample_ids,
      vcf = concat_pbsv_vcf.concatenated_vcf,
      vcf_index = concat_pbsv_vcf.concatenated_vcf_index,
      runtime_attributes = default_runtime_attributes
  }

  call Glnexus.glnexus {
    input:
      cohort_id = family_id + ".joint",
      gvcfs = gvcfs,
      gvcf_indices = gvcf_indices,
      ref_name = ref_map["name"],
      mem_gb = glnexus_mem_gb,
      runtime_attributes = default_runtime_attributes
  }

  call Bcftools.split_vcf_by_sample as split_glnexus {
    input:
      sample_ids = sample_ids,
      vcf = glnexus.vcf,
      vcf_index = glnexus.vcf_index,
      runtime_attributes = default_runtime_attributes
  }

  output {
    Array[File] split_joint_structural_variant_vcfs        = split_pbsv.split_vcfs
    Array[File] split_joint_structural_variant_vcf_indices = split_pbsv.split_vcf_indices
    Array[File] split_joint_small_variant_vcfs             = split_glnexus.split_vcfs
    Array[File] split_joint_small_variant_vcf_indices      = split_glnexus.split_vcf_indices
  }
}
