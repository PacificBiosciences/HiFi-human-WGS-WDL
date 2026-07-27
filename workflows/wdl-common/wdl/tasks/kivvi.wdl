version 1.0

import "../structs.wdl"

task kivvi_kiv2 {
  meta {
    description: "Genotype the LPA KIV2 repeat from HiFi WGS alignments using kivvi."
    outputs: {
      vcf: {
        description: "KIV2 repeat variant VCF"
      },
      vcf_index: {
        description: "Index for KIV2 repeat variant VCF"
      },
      json: {
        description: "KIV2 repeat genotype JSON"
      },
      realigned_bam: {
        description: "KIV2-realigned BAM"
      },
      realigned_bam_index: {
        description: "KIV2-realigned BAM index"
      },
      allele_plot: {
        description: "KIV2 assembled allele plot"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    aligned_bam: {
      description: "Aligned BAM",
      help: "Must be a WGS BAM aligned to GRCh38 at standard depth (20-30x or higher); targeted data is not supported"
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
    out_prefix: {
      description: "Output prefix"
    }
    sensitive: {
      description: "Use a lower fingerprint cutoff for genotyping"
    }
    threads: {
      description: "Number of threads to use"
    }
    mem_gb: {
      description: "Memory to use in GiB"
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
    String out_prefix
    Boolean sensitive = false
    Int threads = 2
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    kivvi --version

    kivvi \
      ~{if sensitive
        then "--sensitive"
        else ""
      } \
      --bam "~{basename(aligned_bam)}" \
      --out . \
      --prefix "~{out_prefix}" \
      --reference "~{basename(ref_fasta)}" \
      kiv2 \
      || echo "kivvi kiv2 failed for ~{basename(aligned_bam)}, presumably due to low coverage." >> messages.txt

    if [ -f "~{out_prefix}.kivvi.kiv2.vcf" ]; then
      bgzip "~{out_prefix}.kivvi.kiv2.vcf"
      tabix --preset vcf "~{out_prefix}.kivvi.kiv2.vcf.gz"
    fi
  >>>

  output {
    File? vcf = "~{out_prefix}.kivvi.kiv2.vcf.gz"
    File? vcf_index = "~{out_prefix}.kivvi.kiv2.vcf.gz.tbi"
    File? json = "~{out_prefix}.kivvi.kiv2.json"
    File? realigned_bam = "~{out_prefix}.kivvi.kiv2.bam"
    File? realigned_bam_index = "~{out_prefix}.kivvi.kiv2.bam.bai"
    File? allele_plot = "~{out_prefix}.kivvi.kiv2.svg"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/kivvi@sha256:55292b4ee7b65645572807a4608832ca6c7e27b0efe2764d671335a7c6d68dca"  # 1.0.0_build2
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task kivvi_d4z4 {
  meta {
    description: "Genotype the D4Z4 repeat from HiFi WGS alignments using kivvi."
    outputs: {
      vcf: {
        description: "D4Z4 repeat variant VCF"
      },
      vcf_index: {
        description: "Index for D4Z4 repeat variant VCF"
      },
      json: {
        description: "D4Z4 repeat genotype JSON"
      },
      realigned_bam: {
        description: "D4Z4-realigned BAM"
      },
      realigned_bam_index: {
        description: "D4Z4-realigned BAM index"
      },
      allele_plot: {
        description: "D4Z4 assembled allele plot"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    aligned_bam: {
      description: "Aligned BAM",
      help: "Must be a WGS BAM aligned to GRCh38 at standard depth (20-30x or higher); targeted data is not supported"
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
    out_prefix: {
      description: "Output prefix"
    }
    sensitive: {
      description: "Use a lower fingerprint cutoff for genotyping"
    }
    threads: {
      description: "Number of threads to use"
    }
    mem_gb: {
      description: "Memory to use in GiB"
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
    String out_prefix
    Boolean sensitive = false
    Int threads = 2
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    kivvi --version

    kivvi \
      ~{if sensitive
        then "--sensitive"
        else ""
      } \
      --bam "~{basename(aligned_bam)}" \
      --out . \
      --prefix "~{out_prefix}" \
      --reference "~{basename(ref_fasta)}" \
      d4z4 \
      || echo "kivvi d4z4 failed for ~{basename(aligned_bam)}, presumably due to low coverage." >> messages.txt

    if [ -f "~{out_prefix}.kivvi.d4z4.vcf" ]; then
      bgzip "~{out_prefix}.kivvi.d4z4.vcf"
      tabix --preset vcf "~{out_prefix}.kivvi.d4z4.vcf.gz"
    fi
  >>>

  output {
    File? vcf = "~{out_prefix}.kivvi.d4z4.vcf.gz"
    File? vcf_index = "~{out_prefix}.kivvi.d4z4.vcf.gz.tbi"
    File? json = "~{out_prefix}.kivvi.d4z4.json"
    File? realigned_bam = "~{out_prefix}.kivvi.d4z4.bam"
    File? realigned_bam_index = "~{out_prefix}.kivvi.d4z4.bam.bai"
    File? allele_plot = "~{out_prefix}.kivvi.d4z4.svg"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/kivvi@sha256:55292b4ee7b65645572807a4608832ca6c7e27b0efe2764d671335a7c6d68dca"  # 1.0.0_build2
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

