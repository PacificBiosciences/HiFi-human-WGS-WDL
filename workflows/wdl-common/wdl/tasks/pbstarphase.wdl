version 1.0

import "../structs.wdl"

task pbstarphase_diplotype {
  meta {
    description: "Run StarPhase to generate diplotype calls and PharmCAT TSV output"
    outputs: {
      summary_json: {
        description: "StarPhase summary"
      },
      pharmcat_tsv: {
        description: "StarPhase PharmCAT TSV"
      }
    }
  }

  parameter_meta {
    out_prefix: {
      description: "Output prefix"
    }
    phased_small_variant_vcf: {
      description: "Phased small variant VCF"
    }
    phased_small_variant_vcf_index: {
      description: "Phased small variant VCF index"
    }
    phased_structural_variant_vcf: {
      description: "Phased structural variant VCF"
    }
    phased_structural_variant_vcf_index: {
      description: "Phased structural variant VCF index"
    }
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Aligned BAM index"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    threads: {
      description: "Number of threads to use"
    }
    mem_gb: {
      description: "Memory to use in GiB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String out_prefix
    File phased_small_variant_vcf
    File phased_small_variant_vcf_index
    File? phased_structural_variant_vcf
    File? phased_structural_variant_vcf_index
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    Int threads = 2
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(phased_small_variant_vcf, "GB") + size(phased_structural_variant_vcf, "GB") + size(aligned_bam, "GB") + size(ref_fasta, "GB") + 50)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{phased_small_variant_vcf}" "~{phased_small_variant_vcf_index}" .
    ~{if defined(phased_structural_variant_vcf)
      then "ln --symbolic --verbose '" + phased_structural_variant_vcf + "' '" + phased_structural_variant_vcf_index + "' ."
      else ""
    }
    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    pbstarphase --version

    pbstarphase diplotype \
      --database /opt/pbstarphase_db.json.gz \
      --reference "~{basename(ref_fasta)}" \
      --vcf "~{basename(phased_small_variant_vcf)}" \
      ~{if defined(phased_structural_variant_vcf)
        then "--sv-vcf '" + basename(select_first([
          phased_structural_variant_vcf
        ])) + "'"
        else ""
      } \
      --bam "~{basename(aligned_bam)}" \
      --output-calls "~{out_prefix}.pbstarphase.json" \
      --pharmcat-tsv "~{out_prefix}.pbstarphase.tsv"
  >>>

  output {
    File summary_json = "~{out_prefix}.pbstarphase.json"
    File pharmcat_tsv = "~{out_prefix}.pbstarphase.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbstarphase@sha256:86aeeb3a9663c22c135d08353de700050e1d268406c1f63832cad5f139cdd390"  # 2.2.0
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

