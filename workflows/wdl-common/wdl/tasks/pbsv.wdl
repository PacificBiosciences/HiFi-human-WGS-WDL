version 1.0

import "../structs.wdl"

task pbsv_discover {
  meta {
    description: "Discover structural variant signatures with `pbsv discover`."
  }

  parameter_meta {
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    trf_bed: {
      name: "Tandem repeat BED used to normalize representation of variation within tandem repeats"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    svsig: {
      name: "Structural variant signature file"
    }
  }

  input {
    File aligned_bam
    File aligned_bam_index

    File? trf_bed

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 10
  Int disk_size = ceil((size(aligned_bam, "GB") + size(trf_bed, "GB")) * 2 + 20)

  String out_prefix = basename(aligned_bam, ".bam")

  command <<<
    set -euo pipefail

    pbsv --version

    pbsv discover \
      --log-level INFO \
      --hifi \
      ~{"--tandem-repeats " + trf_bed} \
      ~{aligned_bam} \
      ~{out_prefix}.svsig.gz
  >>>

  output {
    File svsig = "~{out_prefix}.svsig.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbsv@sha256:3a8529853c1e214809dcdaacac0079de70d0c037b41b43bb8ba7c3fc5f783e26"
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


task pbsv_call {
  meta {
    description: "Call structural variants from signatures with `pbsv call`."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    svsigs: {
      name: "Structural variant signatures"
    }
    sample_count: {
      name: "Number of samples"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    ref_name: {
      name: "Reference name"
    }
    shard_index: {
      name: "Optional shard index"
    }
    regions: {
      name: "Optional shard regions"
    }
    mem_gb: {
      name: "Memory override (GB)"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    vcf: {
      name: "Structural variant VCF"
    }
    vcf_index: {
      name: "Structural variant VCF index"
    }
  }

  input {
    String sample_id
    Array[File] svsigs
    Int? sample_count

    File ref_fasta
    File ref_index
    String ref_name

    Int? shard_index
    Array[String]? regions

    Int mem_gb = if select_first([sample_count, 1]) > 3 then 96 else 64

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int disk_size = ceil((size(svsigs, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  String shard      = if defined(shard_index) then ".~{select_first([shard_index])}" else ""
  String out_prefix = "~{sample_id}.~{ref_name}~{shard}.pbsv"

  command <<<
    set -euo pipefail

    if ~{defined(regions)}; then
    # pbsv has the ability to call SVs by region by using indexed signatures, but
    #   if an svsig.gz file doesn't contain any signatures in the region, then
    #   pbsv crashes. To avoid this, filter the svsig.gz files to only contain
    #   signatures in the regions.
    # This is brittle and likely to break if pbsv discover changes output format.
    # Build a pattern to match; we want headers (e.g., '^#') and signature
    #   records where third column matches the chromosome (e.g., '^.\t.\tchr1\t').
      pattern=$(echo ~{sep=" " select_first([regions, []])} \
        | sed 's/^/^.\\t.\\t/; s/ /\\t\|^.\\t.\\t/g; s/$/\\t/' \
        | echo "^#|""$(</dev/stdin)")

      for svsig in ~{sep=" " svsigs}; do
        filtered_svsig=$(mktemp XXXXXXXXXX.svsig.gz)
        gunzip --stdout "$svsig" \
        | grep --perl-regexp "$pattern" \
        | bgzip --stdout > "${filtered_svsig}" \
        && echo "${filtered_svsig}" >> svsigs.fofn
      done
    else
      cp --verbose ~{write_lines(svsigs)} svsigs.fofn
    fi

    pbsv --version

    pbsv call \
      --hifi \
      --min-sv-length 20 \
      --log-level INFO \
      --num-threads ~{threads} \
      ~{ref_fasta} \
      svsigs.fofn \
      "~{out_prefix}.vcf"

    bgzip --version

    bgzip "~{out_prefix}.vcf"

    tabix --version

    tabix --preset vcf "~{out_prefix}.vcf.gz"
  >>>

  output {
    File vcf       = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbsv@sha256:3a8529853c1e214809dcdaacac0079de70d0c037b41b43bb8ba7c3fc5f783e26"
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
