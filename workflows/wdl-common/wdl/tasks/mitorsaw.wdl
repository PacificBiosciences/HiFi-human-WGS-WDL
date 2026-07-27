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
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    String out_prefix
    Int threads = 4
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

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
    docker: "~{runtime_attributes.container_registry}/mitorsaw@sha256:dc40d68bde4d692254ab0c6959079e61c80b341ef199f71520de0d3c286f72e2"  # 0.2.7_build2
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

