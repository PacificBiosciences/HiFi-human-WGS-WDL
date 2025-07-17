version 1.0

import "../structs.wdl"

task mitorsaw {
  meta {
    description: "Identify and quantify mitochondrial variants and haplotypes from aligned BAM files."
  }

  parameter_meta {
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference index"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attributes"
    }
    vcf: {
      name: "VCF"
    }
    vcf_index: {
      name: "VCF index"
    }
    hap_stats: {
      name: "Haplotype stats"
    }
  }

  input {
    File aligned_bam
    File aligned_bam_index

    File ref_fasta
    File ref_index

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 32
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    mitorsaw haplotype \
      --reference ~{ref_fasta} \
      --bam ~{aligned_bam} \
      --output-vcf ~{out_prefix}.mitorsaw.vcf.gz \
      --output-hap-stats ~{out_prefix}.mitorsaw.json
  >>>

  output {
    File vcf = "~{out_prefix}.mitorsaw.vcf.gz"
    File vcf_index = "~{out_prefix}.mitorsaw.vcf.gz.tbi"
    File hap_stats = "~{out_prefix}.mitorsaw.json"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mitorsaw@sha256:1509dbd7b0a815c7ceb3af52fddc93ef3544ae1858483139450fa0285f8dbe0c"
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
