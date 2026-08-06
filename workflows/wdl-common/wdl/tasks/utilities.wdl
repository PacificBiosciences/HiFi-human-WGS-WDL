version 1.0

import "../structs.wdl"

task split_string {
  meta {
    description: "Split a concatenated string into an array of strings"
    outputs: {
      array: {
        description: "Array of strings obtained by splitting the input string"
      }
    }
  }

  parameter_meta {
    concatenated_string: {
      description: "Concatenated string"
    }
    delimiter: {
      description: "Delimiter"
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
    String concatenated_string
    String delimiter = ","
    Int threads = 1
    Int mem_gb = 1
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 1

  command <<<
    set -euo pipefail
    # Split on delimiter with awk rather than WDL's sub(), whose pattern argument is always a
    # regex: per POSIX/awk semantics, a single-character FS (other than space) is matched
    # literally, so delimiter characters with regex meaning (e.g. "." or "|") still split
    # correctly here. This only holds for single-character delimiters -- a multi-character
    # delimiter would still be interpreted as a regex by awk.
    echo "~{concatenated_string}" | awk -F"~{delimiter}" '{for (i = 1; i <= NF; i++) print $i}'
  >>>

  output {
    #@ except: DeclarationName
    Array[String] array = read_lines(stdout())
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

task consolidate_stats {
  meta {
    description: "Consolidate stats into a TSV file"
    outputs: {
      stats_tsv: {
        description: "TSV file containing consolidated stats"
      },
      messages: {
        description: "Text file containing messages"
      }
    }
  }

  parameter_meta {
    out_prefix: {
      description: "Output prefix"
    }
    stats: {
      description: "Stats"
    }
    msg_array: {
      description: "Array of messages"
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
    String out_prefix
    Array[Array[String]] stats
    #@ except: DeclarationName
    Array[String] msg_array
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 1

  command <<<
    set -euo pipefail
    # flatten the stats map into a TSV format
    # first, like this:
    # stat1 sample1_value1 sample2_value1 ...
    # stat2 sample1_value2 sample2_value2 ...
    # then, transpose the TSV into a TSV with the following format:
    # stat1 stat2 ...
    # sample1_value1 sample1_value2 ...
    # sample2_value1 sample2_value2 ...

    cat << EOF > transpose.py
    import sys
    for outrow in list(zip(*[_.split() for _ in sys.stdin])):
      print('\t'.join(outrow))
    EOF

    python3 ./transpose.py < "~{write_tsv(stats)}" \
    > "~{out_prefix}.stats.txt"

    sed '/^[[:space:]]*$/d' "~{write_lines(msg_array)}" > "~{out_prefix}.messages.txt"
  >>>

  output {
    File stats_tsv = "~{out_prefix}.stats.txt"
    File messages = "~{out_prefix}.messages.txt"
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

