version 1.0

import "../wdl-common/wdl/structs.wdl"
import "assembly_tasks.wdl" as assembly_tasks

workflow assembly {
  meta {
    description: "Given a set of HiFi reads for a human sample, run steps assembly pipeline."
  }

  parameter_meta {
    hifi_reads: {
      name: "HiFi reads (BAMs)"
    }
    sample_id: {
      name: "Sample ID information"
    }
    ref_map_file: {
      name: "TSV containing reference genome information"
    }
  }

  input {
    File hifi_read_bam
    
    String sample_id

    Map[String, String] ref_map

    RuntimeAttributes default_runtime_attributes 

  }
  String pav_sif_path = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/miniwdl_cache/singularity_cache/pav_latest.sif"
  #scatter (hifi_read_bam in hifi_reads) {
  call assembly_tasks.samtools_bam_to_fasta as bam_to_fasta {
    input:
      input_bam          = hifi_read_bam,
      threads            = 32,
      runtime_attributes = default_runtime_attributes
  }
  call assembly_tasks.hifiasm_assembly as hifiasm_assembly {
    input:
      input_fasta        = bam_to_fasta.input_fasta,
      threads            = 32,
      runtime_attributes = default_runtime_attributes
  }
  call assembly_tasks.gfa_to_fa as gfa_to_fa {
    input:
      input_hap1_gfa            = hifiasm_assembly.input_1_asm,
      input_hap2_gfa            = hifiasm_assembly.input_2_asm,
      runtime_attributes = default_runtime_attributes
  }
  call assembly_tasks.pav as pav {
    input:
      ref            = ref_map["hg002_fasta"],
      ref_index            = ref_map["hg002_fasta_index"],
      ref_bgzf_index            = ref_map["hg002_fasta_bgzf_index"],
      input_fa_1     = gfa_to_fa.fasta_hap1,
      input_fa_2     = gfa_to_fa.fasta_hap2,
      sample_name    = sample_id,
      pav_sif_path   = pav_sif_path,
      runtime_attributes = default_runtime_attributes
  }
  call assembly_tasks.liftover as liftover {
    input:
      vcf                = select_first([pav.pav_vcf]),
      chain              = ref_map["hg002_chain"],
      ref                = ref_map["fasta"],
      runtime_attributes = default_runtime_attributes
  }
  call assembly_tasks.filter_liftover as filter_liftover {
    input:
      vcf                = liftover.lifted_vcf,
      runtime_attributes = default_runtime_attributes
  }
 #}

  output {
    #bam to fasta outputs
    File fasta_output = bam_to_fasta.input_fasta
    
    # hifiasm assembly outputs
    File asm_1 = hifiasm_assembly.input_1_asm
    File asm_2 = hifiasm_assembly.input_2_asm

    # pav outputs
    File pav_vcf = select_first([pav.pav_vcf])
    File pav_vcf_index = select_first([pav.pav_vcf_index])

    #liftover outputs
    File lifted_vcf = liftover.lifted_vcf
    File reject_vcf = liftover.reject_vcf 

    #LARGE SV FILTER OUTPUTS
    File large_sv_filtered_vcf = filter_liftover.large_sv_filtered_vcf
    File large_sv_filtered_vcf_index = filter_liftover.large_sv_filtered_vcf_index
  }
}
