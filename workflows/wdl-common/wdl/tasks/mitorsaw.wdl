version 1.0

import "../structs.wdl"

task mitorsaw {
  meta {
    description: "Identify and quantify mitochondrial variants and haplotypes from aligned BAM files"
    outputs: {
      vcf: {
        description: "Mitochondrial variant VCF"
      },
      vcf_index: {
        description: "Index for mitochondrial variant VCF"
      },
      hap_stats: {
        description: "Mitochondrial haplotype statistics"
      }
    }
  }

  parameter_meta {
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
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
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

  Int threads = 4
  Int mem_gb = 32
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{aligned_bam}" .
    ln --symbolic --verbose "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    mitorsaw --version

    mitorsaw haplotype \
      --reference "~{basename(ref_fasta)}" \
      --bam "~{basename(aligned_bam)}" \
      --output-vcf "~{out_prefix}.mitorsaw.vcf.gz" \
      --output-hap-stats "~{out_prefix}.mitorsaw.json"
  >>>

  output {
    File vcf = "~{out_prefix}.mitorsaw.vcf.gz"
    File vcf_index = "~{out_prefix}.mitorsaw.vcf.gz.tbi"
    File hap_stats = "~{out_prefix}.mitorsaw.json"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mitorsaw@sha256:bc1aaa4633d60d753d07b7e079c32767e2283f29b5897be35b59ed97e252d457"  # 0.2.7_build1
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
