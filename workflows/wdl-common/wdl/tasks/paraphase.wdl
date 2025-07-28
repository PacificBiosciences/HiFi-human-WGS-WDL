version 1.0

import "../structs.wdl"

task paraphase {
  meta {
    description: "Haplotype genes in hard to align regions using Paraphase"
  }

  parameter_meta {
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
    sample_id: {
      name: "Sample ID"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    out_json: {
      name: "Paraphase report (JSON)"
    }
    bam: {
      name: "Paraphase realigned BAM"
    }
    bam_index: {
      name: "Paraphase realigned BAM index"
    }
    vcfs_tar: {
      name: "Paraphase VCFs tar"
    }
    msg: {
      name: "Messages"
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

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil(size(aligned_bam, "GB") +size(ref_fasta, "GB") + 20)

  command <<<
    set -eu

    touch messages.txt

    paraphase --version

    paraphase \
      --threads ~{threads} \
      --bam ~{aligned_bam} \
      --reference ~{ref_fasta} \
      --out ./ \
      || echo "Paraphase failed for sample ~{sample_id}.  Check Paraphase logs for details." >> messages.txt

    # tarball the VCFs if they exist
    if ls ~{sample_id}_paraphase_vcfs/*.vcf &> /dev/null; then
      tar --gzip --create --verbose --file ~{sample_id}.paraphase_vcfs.tar.gz ~{sample_id}_paraphase_vcfs/*.vcf
    fi
  >>>

  output {
    File? out_json    = "~{sample_id}.paraphase.json"
    File? bam         = "~{sample_id}.paraphase.bam"
    File? bam_index   = "~{sample_id}.paraphase.bam.bai"
    File? vcfs_tar    = "~{sample_id}.paraphase_vcfs.tar.gz"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/paraphase@sha256:e2f904111a43e8f055681112294e0f05ff2839d9801fc01ac39a17c841016920"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}