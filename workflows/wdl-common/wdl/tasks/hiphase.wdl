version 1.0

import "../structs.wdl"

task hiphase {
  meta {
    description: "Phases VCFs and haplotags BAMs with HiPhase"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    phased_vcf_names: {
      name: "Phased VCF names"
    }
    phased_vcf_index_names: {
      name: "Phased VCF index names"
    }
    aligned_bam: {
      name: "BAM"
    }
    aligned_bam_index: {
      name: "BAM index"
    }
    ref_name: {
      name: "Reference name"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    phased_vcfs: {
      name: "Phased VCFs"
    }
    phased_vcf_indices: {
      name: "Phased VCF indices"
    }
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    phase_stats: {
      name: "Phasing statistics"
    }
    phase_blocks: {
      name: "Phase blocks"
    }
    phase_haplotags: {
      name: "Haplotag information"
    }
    stat_phased_basepairs: {
      name: "Phased basepairs"
    }
    stat_phase_block_ng50: {
      name: "Phase block N50"
    }
    stat_mapped_read_count: {
      name: "Mapped read count"
    }
    stat_mapped_percent: {
      name: "Mapped read percent"
    }
  }

  input {
    String sample_id

    Array[File] vcfs
    Array[File] vcf_indices
    Array[String] phased_vcf_names
    Array[String] phased_vcf_index_names

    File aligned_bam
    File aligned_bam_index

    String ref_name
    File ref_fasta
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = threads * 5
  Int disk_size = ceil(size(vcfs, "GB") + size(ref_fasta, "GB") + size(aligned_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    hiphase --version

    hiphase --threads ~{threads} \
      --sample-name ~{sample_id} \
      ~{sep=" " prefix("--vcf ", vcfs)} \
      ~{sep=" " prefix("--output-vcf ", phased_vcf_names)} \
      --bam ~{aligned_bam} \
      --output-bam ~{sample_id}.~{ref_name}.haplotagged.bam \
      --reference ~{ref_fasta} \
      --summary-file ~{sample_id}.~{ref_name}.hiphase.stats.tsv \
      --blocks-file ~{sample_id}.~{ref_name}.hiphase.blocks.tsv \
      --haplotag-file ~{sample_id}.~{ref_name}.hiphase.haplotags.tsv

    gzip ~{sample_id}.~{ref_name}.hiphase.haplotags.tsv

    # pull the count and percent of mapped reads
    samtools flagstat \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-fmt json \
      ~{sample_id}.~{ref_name}.haplotagged.bam \
      | jq '.["QC-passed reads"]["primary mapped %", "primary mapped"]' \
      > mapping_stats.txt
    sed -n '1p' mapping_stats.txt > mapped_percent.txt
    sed -n '2p' mapping_stats.txt > mapped_read_count.txt

    # pull the phased basepairs and phase block N50
    cat << EOF > get_tsv_stats.py
    import sys, pandas as pd
    df = pd.read_csv('~{sample_id}.~{ref_name}.hiphase.stats.tsv', sep='\t')
    print(df[df['chromosome'] == 'all'][sys.argv[1]].values[0])
    EOF

    python3 get_tsv_stats.py basepairs_per_block_sum > phased_basepairs.txt
    python3 get_tsv_stats.py block_ng50 > phase_block_ng50.txt
  >>>

  output {
    Array  [File] phased_vcfs          = phased_vcf_names
    Array  [File] phased_vcf_indices   = phased_vcf_index_names
    File   haplotagged_bam             = "~{sample_id}.~{ref_name}.haplotagged.bam"
    File   haplotagged_bam_index       = "~{sample_id}.~{ref_name}.haplotagged.bam.bai"
    File   phase_stats                 = "~{sample_id}.~{ref_name}.hiphase.stats.tsv"
    File   phase_blocks                = "~{sample_id}.~{ref_name}.hiphase.blocks.tsv"
    File   phase_haplotags             = "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv.gz"
    String stat_phased_basepairs       = read_string("phased_basepairs.txt")
    String stat_phase_block_ng50       = read_string("phase_block_ng50.txt")
    String stat_mapped_read_count      = read_string("mapped_read_count.txt")
    String stat_mapped_percent         = read_string("mapped_percent.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hiphase@sha256:47fe7d42aea6b1b2e6d3c7401bc35a184464c3f647473d0525c00f3c968b40ad"
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
