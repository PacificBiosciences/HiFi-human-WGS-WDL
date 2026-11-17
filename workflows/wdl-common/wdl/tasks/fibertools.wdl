version 1.0

import "../structs.wdl"

task fibertools_add_nucleosomes {
  meta {
    description: "Add nucleosome calls to a BAM with m6A predictions using fibertools."
    outputs: {
      nucleosome_bam: {
        description: "BAM with nucleosome calls added"
      },
      nucleosome_bam_index: {
        description: "Index for the BAM with nucleosome calls added; only created if the input BAM is aligned"
      }
    }
  }

  parameter_meta {
    m6a_bam: {
      description: "BAM with m6A predictions"
    }
    out_prefix: {
      description: "Output prefix"
    }
    nucleosome_length: {
      description: "Minimum nucleosome length"
    }
    combined_nucleosome_length: {
      description: "Minimum nucleosome length when combining over a single m6A call"
    }
    min_distance_added: {
      description: "Minimum distance needed to add to an already existing nucleosome by crossing an m6A call"
    }
    distance_from_end: {
      description: "Minimum distance from the end of a fiber to call a nucleosome or MSP"
    }
    filter_bit_flag: {
      description: "BAM bit flags to filter on, equivalent to samtools view -F"
    }
    filter_expression: {
      description: "Filtering expression to use for filtering records",
      help: "Supports len() and qual() functions over msp, nuc, m6a, cpg, e.g. \"len(nuc)>150\""
    }
    min_ml_score: {
      description: "Minimum score in the ML tag to use or include in the output"
    }
    uncompressed: {
      description: "Output an uncompressed BAM"
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
    File m6a_bam
    String out_prefix
    Int nucleosome_length = 75
    Int combined_nucleosome_length = 100
    Int min_distance_added = 25
    Int distance_from_end = 45
    Int? filter_bit_flag
    String? filter_expression
    Int min_ml_score = 125
    Boolean uncompressed = false
    Int threads = 8
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(m6a_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{m6a_bam}" .

    ft --version

    ft add-nucleosomes \
      --threads ~{threads} \
      --nucleosome-length ~{nucleosome_length} \
      --combined-nucleosome-length ~{combined_nucleosome_length} \
      --min-distance-added ~{min_distance_added} \
      --distance-from-end ~{distance_from_end} \
      ~{if (defined(filter_bit_flag)) then "--filter '" + filter_bit_flag + "'" else ""} \
      ~{if (defined(filter_expression)) then "--ftx '" + filter_expression + "'" else ""} \
      --ml ~{min_ml_score} \
      ~{true="--uncompressed" false="" uncompressed} \
      "~{basename(m6a_bam)}" \
      "~{out_prefix}.nucleosomes.bam"

    sq_count=$(samtools view -H "~{out_prefix}.nucleosomes.bam" | grep -c '^@SQ' || true)
    if [ "${sq_count}" -gt 0 ]; then
      # shellcheck disable=SC2086
      samtools index \
        ~{if threads > 1 then "--threads '" + (threads - 1) + "'" else ""} \
        "~{out_prefix}.nucleosomes.bam"
    fi
  >>>

  output {
    File nucleosome_bam = "~{out_prefix}.nucleosomes.bam"
    File? nucleosome_bam_index = "~{out_prefix}.nucleosomes.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/fibertools@sha256:bc721a08a447625400efbb152b6be4f0f36cd7b983429f178771fe326fcd3180"  # 0.13.0
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

