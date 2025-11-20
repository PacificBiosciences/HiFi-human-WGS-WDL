version 1.0

import "../structs.wdl"

task pbstarphase_diplotype {
  meta {
    description: "Run PBStarPhase to generate diplotype calls and PharmCAT TSV output."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    phased_small_variant_vcf: {
      name: "Phased small variant VCF file"
    }
    phased_small_variant_vcf_index: {
      name: "Phased small variant VCF index file"
    }
    phased_structural_variant_vcf: {
      name: "Phased structural variant VCF file"
    }
    phased_structural_variant_vcf_index: {
      name: "Phased structural variant VCF index file"
    }
    aligned_bam: {
      name: "Aligned BAM file"
    }
    aligned_bam_index: {
      name: "Aligned BAM index file"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    out_json: {
      name: "PBStarPhase JSON output"
    }
    pharmcat_tsv: {
      name: "PBStarPhase PharmCAT TSV output"
    }
  }

  input {
    String sample_id

    File phased_small_variant_vcf
    File phased_small_variant_vcf_index

    File phased_structural_variant_vcf
    File phased_structural_variant_vcf_index

    File aligned_bam
    File aligned_bam_index

    File ref_fasta
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 16
  Int disk_size = ceil(size(phased_small_variant_vcf, "GB") + size(phased_structural_variant_vcf, "GB") + size(aligned_bam, "GB") + size(ref_fasta, "GB") + 50)

  command <<<
    set -euo pipefail

    pbstarphase --version

    pbstarphase diplotype \
      --database /opt/pbstarphase_db.json.gz \
      --reference ~{ref_fasta} \
      --vcf ~{phased_small_variant_vcf} \
      --sv-vcf ~{phased_structural_variant_vcf} \
      --bam ~{aligned_bam} \
      --output-calls ~{sample_id}.pbstarphase.json \
      --pharmcat-tsv ~{sample_id}.pharmcat.tsv
  >>>

  output {
    File out_json     = "~{sample_id}.pbstarphase.json"
    File pharmcat_tsv = "~{sample_id}.pharmcat.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbstarphase@sha256:b63d7abb718a11c29080cef19026d7c38b8d87df957969dd17b27b15c56de68b"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}