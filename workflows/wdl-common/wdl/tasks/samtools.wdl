version 1.0

import "../structs.wdl"

task samtools_merge {
  meta {
    description: "Merge multiple BAM files into a single BAM file."
  }

  parameter_meta {
    bams: {
      name: "BAMs"
    }
    out_prefix: {
      name: "Output BAM name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_bam: {
      name: "Merged BAM"
    }
    merged_bam_index: {
      name: "Merged BAM index"
    }
  }

  input {
    Array[File] bams

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 4
  Int disk_size = ceil(size(bams, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    samtools merge \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      -o ~{out_prefix}.bam \
      ~{sep=' ' bams}

    samtools index ~{out_prefix}.bam
  >>>

  output {
    File merged_bam       = "~{out_prefix}.bam"
    File merged_bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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

task samtools_fasta {
  meta {
    description: "Convert a BAM file to a FASTA file."
  }

  parameter_meta {
    bam: {
      name: "BAM"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    fasta: {
      name: "FASTA"
    }
  }

  input {
    File bam

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(bam, "GB") * 3.5 + 20)

  String out_prefix = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    samtools --version

    samtools fasta \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{bam} \
    > ~{out_prefix}.fasta
  >>>

  output {
    File fasta = "~{out_prefix}.fasta"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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

task samtools_reset {
  meta {
    description: "Reset SAM/BAM files to remove unwanted tags."
  }

  parameter_meta {
    bam: {
      name: "BAM"
    }
    remove_tags: {
      name: "Tags to remove"
    }
    reject_pg: {
      name: "Reject all PG tags *after* this value"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    out_bam: {
      name: "Output BAM"
    }
  }

  input {
    File bam

    String remove_tags = "HP,PS,PC,SA,mg,mc,mi,rm,fi,fp,ri,rp"
    String reject_pg   = "pbmm2"

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 4
  Int disk_size = ceil(size(bam, "GB") * 3.5 + 20)

  String out_prefix = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    samtools --version

    samtools reset \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --remove-tag ~{remove_tags} \
      --reject-PG ~{reject_pg} \
      -o {out_prefix}.reset.bam \
      ~{bam}
  >>>

  output {
    File out_bam = "~{out_prefix}.reset.bam"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}