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
    phased_vcf: {
      name: "Phased VCF file"
    }
    phased_vcf_index: {
      name: "Phased VCF index file"
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

    File phased_vcf
    File phased_vcf_index

    File aligned_bam
    File aligned_bam_index

    File ref_fasta
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 16
  Int disk_size = ceil(size(phased_vcf, "GB") * 2 + size(aligned_bam, "GB") + size(ref_fasta, "GB") + 50)

  command <<<
    set -euo pipefail

    pbstarphase --version

    pbstarphase diplotype \
      --database /opt/pbstarphase_db.json.gz \
      --reference ~{ref_fasta} \
      --vcf ~{phased_vcf} \
      --bam ~{aligned_bam} \
      --output-calls ~{sample_id}.pbstarphase.json \
      --pharmcat-tsv ~{sample_id}.pharmcat.tsv
  >>>

  output {
    File out_json     = "~{sample_id}.pbstarphase.json"
    File pharmcat_tsv = "~{sample_id}.pharmcat.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbstarphase@sha256:426764fb09eadbc5de8aea2450b5fe55000c1dd7aaa2ea7b0f5438f34ea63e3d"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}