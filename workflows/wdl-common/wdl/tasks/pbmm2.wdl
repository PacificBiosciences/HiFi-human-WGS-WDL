version 1.0

import "../structs.wdl"

task pbmm2_align_wgs {
  meta {
    description: "Align HiFi reads to a reference genome with pbmm2."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    bam: {
      name: "HiFi reads (BAM)"
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
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    bam_stats: {
      name: "BAM stats"
    }
  }

  input {
    String sample_id
    File bam

    File ref_fasta
    File ref_index
    String ref_name

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 24
  Int mem_gb    = ceil(threads * 4)
  Int disk_size = ceil((size(bam, "GB") + size(ref_fasta, "GB")) * 3 + 20)

  String movie = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    pbmm2 --version

    pbmm2 align \
      --num-threads ~{threads} \
      --sort-memory 4G \
      --preset HIFI \
      --sample ~{sample_id} \
      --log-level INFO \
      --sort \
      --unmapped \
      ~{ref_fasta} \
      ~{bam} \
      ~{sample_id}.~{movie}.~{ref_name}.aligned.bam &

    cat << EOF > extract_read_length_and_qual.py
    import math, pysam
    MAX_QV = 60
    save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
    bamin = pysam.AlignmentFile('~{bam}', check_sq=False)
    pysam.set_verbosity(save)  # restore warnings
    for b in bamin:
      errorrate = 1.0 - b.get_tag("rq")
      readqv = MAX_QV if errorrate == 0 else math.floor(-10 * math.log10(errorrate))
      print(f"{b.query_name.split('/')[0]}\t{b.query_name}\t{len(b.query_sequence)}\t{readqv}")
    bamin.close()
    EOF

    python3 ./extract_read_length_and_qual.py \
      | gzip -c > ~{sample_id}.~{movie}.read_length_and_quality.tsv.gz

    wait
  >>>

  output {
    File aligned_bam       = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
    File aligned_bam_index = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
    File bam_stats         = "~{sample_id}.~{movie}.read_length_and_quality.tsv.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:265eef770980d93b849d1ddb4a61ac449f15d96981054e91d29da89943084e0e"
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