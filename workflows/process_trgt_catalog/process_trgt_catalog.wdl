version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/workflows/pbmm2/pbmm2.wdl" as Pbmm2

workflow process_trgt_catalog {
  meta {
    description: "Process a TRGT catalog to identify regions for which to include fail_reads."
    outputs: {
      full_catalog: {
        description: "Repeat catalog for TRGT, flags stripped"
      },
      include_fail_reads_bed: {
        description: "Subset of repeat catalog for which to include fail reads"
      },
      fail_reads_bait_index: {
        description: "MMI index for subset of reference corresponding to regions for which to include fail reads"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    trgt_catalog: {
      description: "Repeat catalog for TRGT"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    File trgt_catalog
    File ref_fasta
    File ref_index
    RuntimeAttributes default_runtime_attributes
  }

  call filter_trgt_catalog { input:
    trgt_catalog = trgt_catalog,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    out_prefix = "fail_reads_subset",
    runtime_attributes = default_runtime_attributes
  }

  if (defined(filter_trgt_catalog.include_fail_reads_bed)) {
    call Pbmm2.create_pbmm2_index { input:
      ref_fasta = select_first([
        filter_trgt_catalog.fasta
      ]),
      runtime_attributes = default_runtime_attributes
    }
  }

  output {
    File full_catalog = filter_trgt_catalog.full_catalog
    File? include_fail_reads_bed = filter_trgt_catalog.include_fail_reads_bed
    File? fail_reads_bait_index = create_pbmm2_index.index
    Array[String] msg = filter_trgt_catalog.msg
  }
}

task filter_trgt_catalog {
  meta {
    description: "Filter and clean a TRGT catalog to identify regions for which to include fail_reads."
    outputs: {
      full_catalog: {
        description: "Repeat catalog for TRGT, flags stripped"
      },
      include_fail_reads_bed: {
        description: "Subset of repeat catalog for which to include fail reads"
      },
      fasta: {
        description: "Output FASTA"
      },
      fasta_index: {
        description: "Output FASTA index"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    trgt_catalog: {
      description: "Repeat catalog for TRGT"
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
    File trgt_catalog
    Int slop_size = 10000
    File ref_fasta
    File ref_index
    String out_prefix
    Int threads = 4
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(trgt_catalog, "GB") + size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .
    ln --symbolic --verbose "~{trgt_catalog}" .

    if gzip --test "~{basename(trgt_catalog)}"; then
      gunzip --stdout "~{basename(trgt_catalog)}" > in.bed
      bed="./in.bed"
    else
      bed="~{basename(trgt_catalog)}"
    fi

    # sanitize general TRGT repeat catalog to remove INCLUDE_FAIL_READS flag
    sed 's/;INCLUDE_FAIL_READS//' "${bed}" > trgt.bed &
    PID=$!

    # Create a catalog with regions that have the INCLUDE_FAIL_READS flag
    # If the flag is not present, remove the file
    grep 'INCLUDE_FAIL_READS' "${bed}" \
    > include_fail_reads.trgt.bed || true

    if [ -s include_fail_reads.trgt.bed ]; then
      echo "INCLUDE_FAIL_READS regions: $(cut -f4 include_fail_reads.trgt.bed | cut -f1 -d';' | cut -f2 -d= | paste -sd, -)" >> messages.txt

      samtools --version
      samtools faidx \
        --region-file <(\
          bedtools slop \
            -b ~{slop_size} \
            -g "~{ref_index}" \
            -i "include_fail_reads.trgt.bed" \
          | awk '{{print $1":"$2"-"$3}}') \
        "~{basename(ref_fasta)}" \
        --write-index \
      > "~{out_prefix}.fasta"
    else
      echo "No repeats in ~{trgt_catalog} contain INCLUDE_FAIL_READS flag. fail_reads will not be used for TRGT." >> messages.txt
      rm --verbose --force include_fail_reads.trgt.bed
    fi

    rm --verbose --force "in.bed"

    wait $PID
  >>>

  output {
    File full_catalog = "trgt.bed"
    File? include_fail_reads_bed = "include_fail_reads.trgt.bed"
    File? fasta = "~{out_prefix}.fasta"
    File? fasta_index = "~{out_prefix}.fasta.fai"
    Array[String] msg = read_lines("messages.txt")
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
    cacheable: true  # !UnknownRuntimeKey
  }
}

