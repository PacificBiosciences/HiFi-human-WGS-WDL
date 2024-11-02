version 1.0

import "../structs.wdl"

task split_string {
  meta {
    description: "Split a concatenated string into an array of strings"
  }

  parameter_meta {
    concatenated_string: {
      name: "Concatenated String"
    }
    delimiter: {
      name: "Delimiter"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    array: {
      name: "Array of strings"
    }
  }

  input {
    String concatenated_string
    String delimiter = ","

    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb  = 1
  Int disk    = 1

  command <<<
    echo '~{sub(concatenated_string, delimiter, "\n")}'
  >>>

  output {
    Array[String] array = read_lines(stdout())
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: "~{disk} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task consolidate_stats {
  meta {
    description: "Consolidate stats into a TSV file"
  }

  parameter_meta {
    id: {
      name: "ID"
    }
    stats: {
      name: "Stats"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    output_tsv: {
      name: "Output TSV"
    }
  }

  input {
    String id

    Map[String, Array[String]] stats

    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb  = 1
  Int disk    = 1

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

    # shellcheck disable=SC2296
    jq -cr 'to_entries[] | [.key, (.value | flatten[])] | @tsv' < ~{write_json(stats)} \
    | python3 ./transpose.py \
    > ~{id}.stats.txt
  >>>

  output {
    File output_tsv = "~{id}.stats.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: "~{disk} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}