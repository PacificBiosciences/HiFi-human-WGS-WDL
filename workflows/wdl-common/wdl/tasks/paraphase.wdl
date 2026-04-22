version 1.0

import "../structs.wdl"

task paraphase {
  meta {
    description: "Haplotype genes in hard to align regions using Paraphase"
    outputs: {
      summary_json: {
        description: "Paraphase summary"
      },
      bam: {
        description: "BAM file of reads realigned by Paraphase"
      },
      bam_index: {
        description: "Index for BAM file of reads realigned by Paraphase"
      },
      vcfs_tar: {
        description: "Paraphase VCFs"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Aligned BAM index"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    sample_id: {
      description: "Sample ID"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  Int threads = 16
  Int mem_gb = 32
  Int disk_size = ceil(size(aligned_bam, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -eu

    touch messages.txt

    ln --symbolic --verbose "~{aligned_bam}" .
    ln --symbolic --verbose "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    paraphase --version

    paraphase \
      --threads ~{threads} \
      --bam "~{basename(aligned_bam)}" \
      --reference "~{basename(ref_fasta)}" \
      --out ./ \
      || echo "Paraphase failed for sample ~{sample_id}.  Check Paraphase logs for details." >> messages.txt

    # tarball the VCFs if they exist
    if ls "~{sample_id}_paraphase_vcfs"/*.vcf &> /dev/null; then
      tar --gzip --create --verbose --file "~{sample_id}.paraphase_vcfs.tar.gz" "~{sample_id}_paraphase_vcfs"/*.vcf
    fi
  >>>

  output {
    File? summary_json = "~{sample_id}.paraphase.json"
    File? bam = "~{sample_id}.paraphase.bam"
    File? bam_index = "~{sample_id}.paraphase.bam.bai"
    File? vcfs_tar = "~{sample_id}.paraphase_vcfs.tar.gz"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/paraphase@sha256:e7e1bd125019211bd1341cf47f884558918c67fe1f3a4da3ea04c26c7aaeebca"  # 3.5.0_build1
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
