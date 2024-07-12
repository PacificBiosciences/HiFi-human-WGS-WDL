version 1.0

import "../structs.wdl"

task pangu_cyp2d6 {
  meta {
    description: "Call CYP2D6 haplotypes using Pangu"
  }

  parameter_meta {
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    mode: {
      name: "Pangu mode",
      default: "capture",
      choices: ["wgs", "amplicon", "capture", "consensus"]
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pangu_json: {
      name: "Pangu JSON"
    }
    pangu_tsv: {
      name: "Pangu TSV"
    }
    fixed_pangu_tsv: {
      name: "Fixed Pangu TSV"
    }
  }

  input {
    File haplotagged_bam
    File haplotagged_bam_index

    String mode

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 12
  Int disk_size = ceil(size(haplotagged_bam, "GB") * 2 + 20)

  String out_prefix = basename(haplotagged_bam, ".bam")

  command <<<
    set -euo pipefail

    pangu \
      -m ~{mode} \
      -p ~{out_prefix}.pangu \
      ~{haplotagged_bam}

    # Fix the pangu output with missing calls for the sample
    awk \
      'BEGIN {{OFS="\t"}} !($2 ~ /\//) {{$2=$2"/[]"}} 1' \
      ~{out_prefix}.pangu_pharmcat.tsv \
      > ~{out_prefix}.pangu_pharmcat_fix.tsv
  >>>

  output {
    File pangu_json      = "~{out_prefix}.pangu_report.json"
    File pangu_tsv       = "~{out_prefix}.pangu_pharmcat.tsv"
    File fixed_pangu_tsv = "~{out_prefix}.pangu_pharmcat_fix.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pangu@sha256:f68320856c709b07f52ca84b1a89494d5146c2660f55fcef591dd3280396ea4a"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
