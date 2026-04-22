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
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String concatenated_string
    String delimiter = ","
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 1
  Int disk_size = 1

  command <<<
    echo '~{sub(concatenated_string, delimiter, "\n")}'
  >>>

  output {
    #@ except: DeclarationName
    Array[String] array = read_lines(stdout())
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
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String out_prefix
    Array[Array[String]] stats
    #@ except: DeclarationName
    Array[String] msg_array
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 1
  Int disk_size = 1

  command <<<
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
