version 1.0

import "../structs.wdl"

task owl {
  meta {
    description: "Microsatellite instability (MSI) analysis for HiFi data"
    outputs: {
      profile: {
        description: "Intermediate genotyping output"
      },
      score: {
        description: "Sample summary including MSI score"
      }
    }
  }

  parameter_meta {
    out_prefix: {
      description: "Owl output prefix"
    }
    haplotagged_bam: {
      description: "Aligned, haplotype-tagged BAM"
    }
    haplotagged_bam_index: {
      description: "Aligned, haplotype-tagged BAM index"
    }
    marker_bed: {
      description: "BED of microsatellite marker regions, motifs"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String out_prefix
    File haplotagged_bam
    File haplotagged_bam_index
    File marker_bed
    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb = 16
  Int disk_size = floor(size(haplotagged_bam, "GB")) + 20

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{haplotagged_bam}" .
    ln --symbolic --verbose "~{haplotagged_bam_index}" .
    ln --symbolic --verbose "~{marker_bed}" .

    owl --version
    owl profile \
      --regions "~{basename(marker_bed)}" \
      --bam "~{basename(haplotagged_bam)}" \
      > "~{out_prefix}.owl-profile.txt"

    owl score \
      --file "~{out_prefix}.owl-profile.txt" \
      --prefix "~{out_prefix}"
  >>>

  output {
    File profile = "~{out_prefix}.owl-profile.txt"
    File score = "~{out_prefix}.owl-scores.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/owl@sha256:341b9610a058786da55c208fe6469793857207c3dee58ee2b734bba5f1e4b711"  # 0.4.0_build1
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
