version 1.0

import "../structs.wdl"

task trgt {
  meta {
    description: "Genotype tandem repeats from aligned reads using TRGT."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    trgt_bed: {
      name: "TRGT tandem repeat catalog BED"
    }
    out_prefix: {
      name: "Output prefix"
    }
    max_depth: {
      name: "Maximum locus depth"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    bam: {
      name: "TRGT spanning reads BAM"
    }
    bam_index: {
      name: "TRGT spanning reads BAM index"
    }
    vcf: {
      name: "TRGT repeats VCF"
    }
    vcf_index: {
      name: "TRGT repeats VCF index"
    }
  }

  input {
    String sample_id
    String? sex

    File aligned_bam
    File aligned_bam_index

    File ref_fasta
    File ref_index

    File trgt_bed

    String out_prefix

    Int max_depth = 50

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 32
  Int mem_gb    = 32
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  Int samtools_sort_threads = 8

  String karyotype = if select_first([sex, "FEMALE"]) == "MALE" then "XY" else "XX"

  command <<<
    set -euo pipefail

    echo ~{if defined(sex) then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for TRGT."}

    trgt --version

    trgt genotype \
      --threads ~{threads} \
      --karyotype ~{karyotype} \
      --genome ~{ref_fasta} \
      --repeats ~{trgt_bed} \
      --reads ~{aligned_bam} \
      --max-depth ~{max_depth} \
      --output-prefix ~{out_prefix}.trgt

    bcftools --version

    bcftools sort \
      --output-type z \
      --output ~{out_prefix}.trgt.sorted.vcf.gz \
      ~{out_prefix}.trgt.vcf.gz

    bcftools index \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --tbi \
      ~{out_prefix}.trgt.sorted.vcf.gz

    samtools --version

    samtools sort \
      --threads ~{samtools_sort_threads} \
      -m 400M \
      -o ~{out_prefix}.trgt.spanning.sorted.bam \
      ~{out_prefix}.trgt.spanning.bam

    samtools index \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{out_prefix}.trgt.spanning.sorted.bam

    bcftools view --no-header --exclude-uncalled \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{out_prefix}.trgt.sorted.vcf.gz \
      | wc --lines > genotyped_count.txt || echo "0" > genotyped_count.txt

    bcftools view --no-header --uncalled \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{out_prefix}.trgt.sorted.vcf.gz \
    | wc --lines > uncalled_count.txt || echo "0" > uncalled_count.txt
  >>>

  output {
    File   bam                  = "~{out_prefix}.trgt.spanning.sorted.bam"
    File   bam_index            = "~{out_prefix}.trgt.spanning.sorted.bam.bai"
    File   vcf                  = "~{out_prefix}.trgt.sorted.vcf.gz"
    File   vcf_index            = "~{out_prefix}.trgt.sorted.vcf.gz.tbi"
    String stat_genotyped_count = read_string("genotyped_count.txt")
    String stat_uncalled_count  = read_string("uncalled_count.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:be7e6ef589a31f4de5d2ed4725dfb34b4b23cb9a440577b606e8f7bfee06526b"
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

task trgt_merge {
  meta {
    description: "Merge results multiple samples analyzed with TRGT."
  }

  parameter_meta {
    vcfs: {
      name: "TRGT VCFs"
    }
    vcf_indices: {
      name: "TRGT VCF indices"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_vcf: {
      name: "Merged TRGT repeats VCF"
    }
    merged_vcf_index: {
      name: "Merged TRGT repeats VCF index"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices

    File ref_fasta
    File ref_index

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = ceil((size(vcfs, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    trgt --version

    trgt merge \
      --vcf ~{sep=" " vcfs} \
      --genome ~{ref_fasta} \
      --output ~{out_prefix}.vcf.gz

    bcftools index \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --tbi \
      ~{out_prefix}.vcf.gz
  >>>

  output {
    File merged_vcf       = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:be7e6ef589a31f4de5d2ed4725dfb34b4b23cb9a440577b606e8f7bfee06526b"
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

task coverage_dropouts {
  meta {
    description: "Get coverage dropouts from aligned reads."
  }

  parameter_meta {
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    trgt_bed: {
      name: "TRGT tandem repeat catalog BED"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    dropouts: {
      name: "TRGT regions with coverage dropouts"
    }
  }

  input {
    File aligned_bam
    File aligned_bam_index

    File trgt_bed

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil((size(aligned_bam, "GB")) + 20)

  command <<<
    set -euo pipefail

    # Get coverage dropouts
    check_trgt_coverage.py \
      ~{trgt_bed} \
      ~{aligned_bam} \
    > ~{out_prefix}.trgt.dropouts.txt
  >>>

  output {
    File dropouts = "~{out_prefix}.trgt.dropouts.txt"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:be7e6ef589a31f4de5d2ed4725dfb34b4b23cb9a440577b606e8f7bfee06526b"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}