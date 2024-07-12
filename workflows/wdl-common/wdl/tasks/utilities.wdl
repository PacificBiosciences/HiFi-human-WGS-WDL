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

  Int threads = 2
  Int mem_gb  = 2
  Int disk    = 4

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

task read_pbsv_splits {
  meta {
    description: "Read a JSON file containing an array of arrays of strings; equivalent to read_json() but allows String->File coercion"
  }

  parameter_meta {
    pbsv_splits_file: {
      name: "PBSV splits file"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    splits: {
      name: "Array of arrays of strings"
    }
  }
  
  input {
    File pbsv_splits_file

    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb  = 2
  Int disk    = 4

  command <<<
    cat ~{pbsv_splits_file} > out.json
  >>>

  output {
    Array[Array[String]] splits = read_json("out.json")
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