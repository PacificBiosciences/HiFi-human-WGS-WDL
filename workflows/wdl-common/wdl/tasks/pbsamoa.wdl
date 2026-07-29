version 1.0

import "../structs.wdl"

task pbsamoa_merge {
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
      description: "BAMs to merge, must be coordinate-sorted"
    }
    out_prefix: {
      description: "Output BAM prefix"
    }
    compression: {
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
    Int compression = 6
    Int threads = 32
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(bams, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    # pbsamoa auto-detects coordinate- vs queryname- vs unsorted merge strategy from each
    # input BAM's @HD SO: tag; inputs here are always coordinate-sorted, which takes the
    # multithreaded heap-merge path below (an all-unsorted/unknown input set would instead
    # silently fall back to a single-threaded concat path that ignores --compress-threads/
    # --decode-threads entirely).
    #shellcheck disable=SC2086
    pbsamoa merge \
      --compress-threads ~{ceil(threads * 3 / 4)} \
      --decode-threads ~{floor(threads * 1 / 4)} \
      --memory ~{mem_gb}G \
      --compression ~{compression} \
      --bai \
      "~{out_prefix}.bam" \
      ~{sep=" " bams}
  >>>

  output {
    File merged_bam = "~{out_prefix}.bam"
    File merged_bam_index = "~{out_prefix}.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbsamoa@sha256:381891341e4d33ea9b3ca1c8a7cac7d86d4feea282bc8d12b6201a77cbd4216e"  # 20250702_build1
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

