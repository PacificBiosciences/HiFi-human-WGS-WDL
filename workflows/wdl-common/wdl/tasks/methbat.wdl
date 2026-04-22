version 1.0

import "../structs.wdl"

task methbat {
  meta {
    description: "Profile methylation regions with MethBat"
    outputs: {
      profile: {
        description: "MethBat 5mCpG profile"
      },
      stat_methbat_methylated_count: {
        description: "Number of profiled regions labeled as methylated"
      },
      stat_methbat_unmethylated_count: {
        description: "Number of profiled regions labeled as unmethylated"
      },
      stat_methbat_asm_count: {
        description: "Number of profiled regions labeled as having allele-specific methylation"
      }
    }
  }

  parameter_meta {
    sample_prefix: {
      description: "Sample prefix"
    }
    methylation_pileup_beds: {
      description: "Methylation pileup BED files"
    }
    region_tsv: {
      description: "Regions TSV"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_prefix
    Array[File] methylation_pileup_beds
    File region_tsv
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb = 4
  Int disk_size = ceil((size(methylation_pileup_beds, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    for i in ~{sep=" " methylation_pileup_beds}; do
      ln --symbolic --verbose "${i}" .
    done

    methbat --version

    methbat profile \
      --input-prefix "~{sample_prefix}" \
      --input-regions "~{region_tsv}" \
      --output-region-profile "~{out_prefix}.methbat.profile.tsv"

    # count for three most interesting methylation summary_labels
    awk '$5=="Methylated" {print}' "~{out_prefix}.methbat.profile.tsv" \
    | wc -l > methylated_count.txt
    awk '$5=="Unmethylated" {print}' "~{out_prefix}.methbat.profile.tsv" \
    | wc -l > unmethylated_count.txt
    awk '$5=="AlleleSpecificMethylation" {print}' "~{out_prefix}.methbat.profile.tsv" \
    | wc -l > asm_count.txt
  >>>

  output {
    File profile = "~{out_prefix}.methbat.profile.tsv"
    String stat_methbat_methylated_count = read_string("methylated_count.txt")
    String stat_methbat_unmethylated_count = read_string("unmethylated_count.txt")
    String stat_methbat_asm_count = read_string("asm_count.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:0a5363af6a8ba7670e87fa7d6754bc3afe85dfb73e5deabc4c21913323a5bdd9"  # 0.17.0_build1
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
