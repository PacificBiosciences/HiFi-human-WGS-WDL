version 1.0

import "../../structs.wdl"

workflow get_pbsv_splits {
  meta {
    description: "Parse the pbsv_splits file and return the splits"
  }

  parameter_meta {
    pbsv_splits_file: {
      name: "File describing how pbsv should be parallelized"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pbsv_splits: {
      name: "Array of arrays of strings, each outer array is a list of chromosomes to be processed together by pbsv"
    }
  }

  input {
    File pbsv_splits_file

    RuntimeAttributes default_runtime_attributes
  }

  if (default_runtime_attributes.backend == "AWS-HealthOmics") {
    call read_pbsv_splits {
      input:
        pbsv_splits_file   = pbsv_splits_file,
        runtime_attributes = default_runtime_attributes
    }
  }

  if (default_runtime_attributes.backend != "AWS-HealthOmics") {
    Array[Array[String]] native_splits = read_json(pbsv_splits_file)
  }

  output {
    Array[Array[String]] pbsv_splits = select_first([read_pbsv_splits.splits, native_splits])
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

  Int threads = 1
  Int mem_gb  = 1
  Int disk    = 1

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