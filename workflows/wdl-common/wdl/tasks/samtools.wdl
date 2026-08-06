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
    compression_level: {
      description: "Compression level for the output BAM"
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
    Array[File] bams
    String out_prefix
    Int compression_level = -1
    Int threads = 16
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

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
      -l ~{compression_level} \
      -c -p \
      --write-index \
      -o "~{out_prefix}.bam##idx##~{out_prefix}.bam.bai" \
      ~{sep=" " bams}
  >>>

  output {
    File merged_bam = "~{out_prefix}.bam"
    File merged_bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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
    exclude_mask: {
      description: "Exclude reads with this optional mask, default is to include all reads"
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
    File bam
    Int? exclude_mask
    Int threads = 16
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

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
      ~{if defined(exclude_mask)
        then "--exclude-flags '" + exclude_mask + "'"
        else ""
      } \
      "~{basename(bam)}" \
    > "~{out_prefix}.fasta"
  >>>

  output {
    File fasta = "~{out_prefix}.fasta"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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
    File bed
    Int slop_size = 10000
    File ref_fasta
    File ref_index
    String out_prefix
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bed}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    samtools --version

    samtools faidx \
      --region-file <(\
        bedtools slop \
          -b ~{slop_size} \
          -g "~{ref_index}" \
          -i "~{basename(bed)}" \
        | awk '{{print $1":"$2"-"$3}}') \
      "~{basename(ref_fasta)}" \
      --write-index \
    > "~{out_prefix}.fasta"
    samtools faidx "~{out_prefix}.fasta"
  >>>

  output {
    File fasta = "~{out_prefix}.fasta"
    File fasta_index = "~{out_prefix}.fasta.fai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

task subset_bam {
  meta {
    description: "Create a subset of an aligned BAM file based on a BED file"
    outputs: {
      bam: {
        description: "Output BAM"
      },
      bam_index: {
        description: "Output BAM index"
      },
      count_passed: {
        description: "Count of passed records"
      },
      count_failed: {
        description: "Count of failed records"
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
    File bed
    File aligned_bam
    File aligned_bam_index
    File ref_index
    Int slop_size = 10000
    String out_prefix
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bed}" .
    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
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
      --save-counts "~{out_prefix}.counts.json" \
        --write-index \
      --output "~{out_prefix}.bam##idx##~{out_prefix}.bam.bai" \
      "~{basename(aligned_bam)}"

    jq '.records_filter_accepted' "~{out_prefix}.counts.json" > "~{out_prefix}.count_passed.txt"
    jq '.records_filter_rejected' "~{out_prefix}.counts.json" > "~{out_prefix}.count_failed.txt"
  >>>

  output {
    File bam = "~{out_prefix}.bam"
    File bam_index = "~{out_prefix}.bam.bai"
    String count_passed = read_string("~{out_prefix}.count_passed.txt")
    String count_failed = read_string("~{out_prefix}.count_failed.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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
    File bam
    String remove_tags = "fi,ri,fp,rp,ip,pw,HP,PS,PC,mc,mg,mi,rm"
    String reject_pg = "pbmm2"
    Int threads = 16
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

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
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

