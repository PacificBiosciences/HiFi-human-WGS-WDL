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
    genome: {
      description: "Genome reference build",
      choices: [
        "38",
        "37",
        "19",
        "chm13"
      ]
    }
    sample_id: {
      description: "Sample ID"
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
    String genome = "38"
    String sample_id
    Int threads = 8
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  # Paraphase's region list (~162 segmental-duplication regions, per
  # https://github.com/PacificBiosciences/paraphase/blob/main/docs/regions.md)
  # and its default depth-normalization behavior (it assumes uniform depth
  # across the genome unless run with --targeted, which this task does not
  # use) are keyed off --genome. Only the "38" default has test coverage here;
  # the --config region catalog (default `data/38/config.yaml` inside the
  # container, per `paraphase --help`) still isn't exposed as an input, so
  # non-38 builds beyond Paraphase's own bundled configs aren't supported yet.
  Int disk_size = ceil(size(aligned_bam, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    paraphase --version

    # Paraphase can catch a panic for an individual region internally and
    # continue on to the rest, exiting 0 overall (observed: 2 of 164 regions
    # panicking with "index out of bounds" at a very high thread count, well
    # beyond what --threads is ever set to here -- see TEST_RECORD.md Known
    # issues). The `||` fallback below only catches a total command failure,
    # so those per-region failures need to be pulled out of its own log
    # separately.
    paraphase \
      --threads ~{threads} \
      --bam "~{basename(aligned_bam)}" \
      --reference "~{basename(ref_fasta)}" \
      --genome "~{genome}" \
      --out ./ \
      2>&1 | tee paraphase.log \
      || echo "Paraphase failed for sample ~{sample_id}.  Check Paraphase logs for details." >> messages.txt

    grep "Gene processing failed for" paraphase.log >> messages.txt || true

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
    docker: "~{runtime_attributes.container_registry}/paraphase@sha256:665ac4fefcef92e0023395eae9e0ca3f6e65d1228c49afba3ad3f8e7b3f3d1cb"  # 4.0.0_build2
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

