version 1.0

import "../structs.wdl"

task parabricks_pbmm2 {
  meta {
    description: "Align HiFi reads to a reference genome using NVIDIA Parabricks pbmm2."
    outputs: {
      aligned_bam: {
        description: "Aligned BAM"
      },
      aligned_bam_index: {
        description: "Aligned BAM index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    bam: {
      description: "Unaligned BAM"
    }
    pbmm2_index: {
      description: "pbmm2 reference index"
    }
    ref_name: {
      description: "Reference name"
    }
    parabricks_version: {
      description: "Parabricks version"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    File pbmm2_index
    String ref_name
    String parabricks_version = "4.7.0-1"
    RuntimeAttributes runtime_attributes
  }

  String docker_image = if (runtime_attributes.backend == "AWS-OMICS")
    then runtime_attributes.container_registry
    else "nvcr.io/nvidia/clara" + "/clara-parabricks:~{parabricks_version}"

  Int threads = 128
  Int mem_gb = 500
  Int gpuCount = 2

  Int max_queue_reads = 5000000
  Int max_queue_chunks = 10000

  Int disk_size = ceil(size(bam, "GB") * 2 + size(pbmm2_index, "GB") + 70)

  String movie = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    # shellcheck disable=SC2034
    TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456

    ln --symbolic --verbose "~{bam}" .
    ln --symbolic --verbose "~{pbmm2_index}" .

    /usr/local/parabricks/pbrun minimap2 \
      --verbose \
      --pbmm2 --pbmm2-unmapped --eqx \
      --num-threads ~{threads} \
      --gpusort --gpuwrite \
      --num-gpus ~{gpuCount} \
      --max-queue-reads ~{max_queue_reads} \
      --max-queue-chunks ~{max_queue_chunks} \
      --process-large-alignments-on-gpu \
      --ref "~{basename(pbmm2_index)}" \
      --in-bam "~{basename(bam)}" \
      --out-bam "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
  >>>

  output {
    File aligned_bam = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
    File aligned_bam_index = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    gpuCount: gpuCount
    gpuType: runtime_attributes.gpuType
    acceleratorCount: gpuCount  # !UnknownRuntimeKey
    acceleratorType: runtime_attributes.gpuType  # !UnknownRuntimeKey
    nvidiaDriverVersion: "535.230.02"  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform  # !UnknownRuntimeKey
  }
}
