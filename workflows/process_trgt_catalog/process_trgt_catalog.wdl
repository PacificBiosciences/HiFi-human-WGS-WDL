version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/samtools.wdl" as Samtools

workflow process_trgt_catalog {
  meta {
    description: "Process a TRGT catalog to identify regions for which to include fail_reads."
  }

  parameter_meta {
    trgt_catalog: {
      name: "Repeat catalog for TRGT; BED"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    default_runtime_attributes: {
      name: "Default runtime attribute structure"
    }
    full_catalog: {
      name: "Repeat catalog for TRGT, flags stripped; BED"
    }
    include_fail_reads_bed: {
      name: "Subset of repeat catalog for which to include fail reads; BED"
    }
    fail_reads_bait_fasta: {
      name: "FASTA of reference sequences for baiting fail reads; FASTA"
    }
    fail_reads_bait_index: {
      name: "Index of reference sequences for baiting fail reads; FASTA index"
    }
    msg: {
      name: "Array of messages"
    }
  }

  input {
    File trgt_catalog

    File ref_fasta
    File ref_index

    RuntimeAttributes default_runtime_attributes
  }

  call filter_trgt_catalog {
    input:
      trgt_catalog = trgt_catalog,
      runtime_attributes = default_runtime_attributes
  }

  if (defined(filter_trgt_catalog.include_fail_reads_bed)) {
    call Samtools.subset_reference {
      input:
        bed = select_first([filter_trgt_catalog.include_fail_reads_bed]),
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        out_prefix = "fail_reads_subset",
        runtime_attributes = default_runtime_attributes
    }
  }

  output {
    File  full_catalog           = filter_trgt_catalog.full_catalog
    File? include_fail_reads_bed = filter_trgt_catalog.include_fail_reads_bed
    File? fail_reads_bait_fasta  = subset_reference.fasta
    File? fail_reads_bait_index  = subset_reference.fasta_index
    Array[String] msg            = filter_trgt_catalog.msg
  }

}

task filter_trgt_catalog {
  meta {
    description: "Filter and clean a TRGT catalog to identify regions for which to include fail_reads."
  }

  parameter_meta {
    trgt_catalog: {
      name: "Repeat catalog for TRGT; BED"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    full_catalog: {
      name: "Repeat catalog for TRGT, flags stripped; BED"
    }
    include_fail_reads_bed: {
      name: "Subset of repeat catalog for which to include fail reads; BED"
    }
    msg: {
      name: "Array of messages"
    }
  }

  input {
    File trgt_catalog

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(trgt_catalog, "GB") + 20)

  command <<<
    set -eu

    touch messages.txt

    if gzip --test ~{trgt_catalog}; then
      gunzip --stdout ~{trgt_catalog} > in.bed
      bed="./in.bed"
    else
      bed="~{trgt_catalog}"
    fi

    # sanitize general TRGT repeat catalog to remove INCLUDE_FAIL_READS flag
    sed 's/;INCLUDE_FAIL_READS//' ${bed} > trgt.bed

    # Create a catalog with regions that have the INCLUDE_FAIL_READS flag and sanitize to remove INCLUDE_FAIL_READS flag
    # If the flag is not present, remove the file
    grep 'INCLUDE_FAIL_READS' ${bed} \
    | sed 's/;INCLUDE_FAIL_READS//' \
    > include_fail_reads.trgt.bed || true

    if [ ! -s include_fail_reads.trgt.bed ]; then
      echo "No repeats in ~{trgt_catalog} contain INCLUDE_FAIL_READS flag. fail_reads will not be used for TRGT." >> messages.txt
      rm include_fail_reads.trgt.bed
    fi
  >>>

  output {
    File  full_catalog           = "trgt.bed"
    File? include_fail_reads_bed = "include_fail_reads.trgt.bed"
    Array[String] msg            = read_lines("messages.txt")
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
