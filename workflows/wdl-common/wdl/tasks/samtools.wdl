version 1.0

import "../structs.wdl"

task samtools_merge {
  meta {
    description: "Merge multiple BAM files into a single BAM file"
    outputs: {
      merged_bam: {
        description: "Merged BAM"
      },
      merged_bam_index: {
        description: "Merged BAM index"
      }
    }
  }

  parameter_meta {
    bams: {
      description: "BAMs to merge"
    }
    out_prefix: {
      description: "Output BAM prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] bams
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(bams, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    samtools --version

    # shellcheck disable=SC2086
    samtools merge \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      -c -p \
      -o "~{out_prefix}.bam" \
      ~{sep=" " bams}

    samtools index "~{out_prefix}.bam"
  >>>

  output {
    File merged_bam = "~{out_prefix}.bam"
    File merged_bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
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

task samtools_fasta {
  meta {
    description: "Convert a BAM file to a FASTA file"
    outputs: {
      fasta: {
        description: "FASTA"
      }
    }
  }

  parameter_meta {
    bam: {
      description: "BAM"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File bam
    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(bam, "GB") * 3.5 + 20)

  String out_prefix = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bam}" .

    samtools --version

    samtools fasta \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      "~{basename(bam)}" \
    > "~{out_prefix}.fasta"
  >>>

  output {
    File fasta = "~{out_prefix}.fasta"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
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

task subset_reference {
  meta {
    description: "Create a subset of a reference FASTA file based on a BED file"
    outputs: {
      fasta: {
        description: "Output FASTA"
      },
      fasta_index: {
        description: "Output FASTA index"
      }
    }
  }

  parameter_meta {
    bed: {
      description: "Regions BED"
    }
    slop_size: {
      description: "Size in base pairs to extend the regions in the BED file"
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
    File bed
    Int slop_size = 10000
    File ref_fasta
    File ref_index
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bed}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    samtools --version

    samtools faidx \
      --region-file <(\
        bedtools slop \
          -b ~{slop_size} \
          -g "~{ref_index}" \
          -i "~{basename(bed)}" \
        | awk '{{print $1":"$2"-"$3}}') \
      "~{basename(ref_fasta)}" \
    > "~{out_prefix}.fasta"
    samtools faidx "~{out_prefix}.fasta"
  >>>

  output {
    File fasta = "~{out_prefix}.fasta"
    File fasta_index = "~{out_prefix}.fasta.fai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task subset_bam {
  meta {
    description: "Create a subset of an aligned BAM file based on a BED file"
    outputs: {
      bam: {
        description: "Output BAM"
      },
      bam_index: {
        description: "Output BAM index"
      }
    }
  }

  parameter_meta {
    bed: {
      description: "Regions BED"
    }
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Aligned BAM index"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    slop_size: {
      description: "Size in base pairs to extend the regions in the BED file"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
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

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bed}" .
    ln --symbolic --verbose "~{aligned_bam}" .
    ln --symbolic --verbose "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_index}" .

    samtools --version

    samtools view \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --targets-file <(\
        bedtools slop \
          -b ~{slop_size} \
          -g "~{basename(ref_index)}" \
          -i "~{basename(bed)}") \
      --output "~{out_prefix}.bam" \
      "~{basename(aligned_bam)}"

    samtools index \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } "~{out_prefix}.bam"
  >>>

  output {
    File bam = "~{out_prefix}.bam"
    File bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task samtools_reset {
  meta {
    description: "Reset SAM/BAM files to remove unwanted tags"
    outputs: {
      reset_bam: {
        description: "Output BAM"
      }
    }
  }

  parameter_meta {
    bam: {
      description: "BAM"
    }
    remove_tags: {
      description: "Tags to remove"
    }
    reject_pg: {
      description: "Reject all PG tags *after* this value"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File bam
    String remove_tags = "HP,PS,PC,SA,mg,mc,mi,rm,fi,fp,ri,rp"
    String reject_pg = "pbmm2"
    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(bam, "GB") * 3.5 + 20)

  String out_prefix = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bam}" .

    samtools --version

    samtools reset \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --remove-tag "~{remove_tags}" \
      --reject-PG "~{reject_pg}" \
      -o "~{out_prefix}.reset.bam" \
      "~{basename(bam)}"
  >>>

  output {
    File reset_bam = "~{out_prefix}.reset.bam"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
