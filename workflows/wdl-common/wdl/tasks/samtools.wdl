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
  Int mem_gb    = 16
  Int disk_size = ceil(size(bams, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    samtools merge \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      -c -p \
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
  Int mem_gb    = 16
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

task subset_reference {
  meta {
    description: "Create a subset of a reference FASTA file based on a BED file."
  }

  parameter_meta {
    bed: {
      name: "Regions BED"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    fasta: {
      name: "Output BAM"
    }
    fasta_index: {
      name: "Output BAM index"
    }
  }

  input {
    File bed

    Int slop_size = 10000

    File ref_fasta
    File ref_index

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 4
  Int disk_size = ceil(size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    samtools faidx \
      --region-file <(\
        bedtools slop \
          -b ~{slop_size} \
          -g ~{ref_index} \
          -i ~{bed} \
        | awk '{{print $1":"$2"-"$3}}') \
      ~{ref_fasta} \
    > ~{out_prefix}.fasta

    samtools faidx ~{out_prefix}.fasta
  >>>

  output {
    File fasta = "~{out_prefix}.fasta"
    File fasta_index = "~{out_prefix}.fasta.fai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task subset_bam {
  meta {
    description: "Create a subset of an aligned BAM file based on a BED file."
  }

  parameter_meta {
    bed: {
      name: "Regions BED"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    slop_size: {
      name: "Size in base pairs to extend the regions in the BED file."
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    bam: {
      name: "Output BAM"
    }
    bam_index: {
      name: "Output BAM index"
    }
  }

  input {
    File bed

    File aligned_bam
    File aligned_bam_index

    File ref_index

    Int slop_size = 10000

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 4
  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    samtools view \
      --threads ~{if threads > 1 then threads - 1 else 1} \
      --targets-file <(\
        bedtools slop \
          -b ~{slop_size} \
          -g ~{ref_index} \
          -i ~{bed}) \
      --output ~{out_prefix}.bam \
      ~{aligned_bam}

    samtools index --threads ~{if threads > 1 then threads - 1 else 1} ~{out_prefix}.bam
  >>>

  output {
    File bam = "~{out_prefix}.bam"
    File bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
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
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}