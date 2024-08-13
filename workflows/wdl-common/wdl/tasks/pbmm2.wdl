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
    strip_kinetics: {
      name: "Strip kinetics tags"
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

    Boolean strip_kinetics = true

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 24
  Int mem_gb    = ceil(threads * 4)
  Int disk_size = ceil(size(bam, "GB") * 3 + size(ref_fasta, "GB") + 70)

  String movie = basename(bam, ".bam")

  command <<<
    set -euo pipefail

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
      | gzip -c > ~{sample_id}.~{movie}.read_length_and_quality.tsv.gz &
    BAM_STATS_PID=$!

    samtools --version

    ALIGNED=0

    # is the BAM already aligned?
    if (samtools view -H ~{bam} | grep "^@PG" | grep -qE "ID:(pbmm2|minimap2|ngmlr)"); then
      echo "Input ~{basename(bam)} is already aligned.  Alignments and and haplotype tags will be stripped."
      ALIGNED=1
    fi

    # does the BAM contain basemod tags?
    if ! (samtools view ~{bam} | head -n 50 | cut -f12- | grep -qE "MM:|Mm:|ML:|Ml:"); then
      echo "Input ~{basename(bam)} does not contain base modification tags.  5mCpG pileups will not be generated."
    fi

    # does the BAM contain kinetics tags?
    if (samtools view ~{bam} | head -n 50 | cut -f12- | grep -qE "fi:|fp:|ri:|rp:"); then
      echo "Input ~{basename(bam)} contains consensus kinetics tags."
      if [ "~{strip_kinetics}" = true ]; then
        echo "Kinetics will be stripped from the output."
      fi
    fi

    pbmm2 --version

    pbmm2 align \
      --num-threads ~{threads} \
      --sort-memory 4G \
      --preset HIFI \
      --sample ~{sample_id} \
      --log-level INFO \
      --sort \
      ~{true='--strip' false='' strip_kinetics} \
      --unmapped \
      ~{ref_fasta} \
      ~{bam} \
      aligned.bam

    if [ "$ALIGNED" -eq 1 ]; then
      # remove haplotype tags
      samtools view \
        ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
        --bam --no-PG \
        --remove-tag HP,PS,PC,SA \
        -o ~{sample_id}.~{movie}.~{ref_name}.aligned.bam \
        aligned.bam &&
          \ rm aligned.bam
      samtools index \
        ~{if threads > 1 then "-@ " + (threads - 1) else ""} \
        ~{sample_id}.~{movie}.~{ref_name}.aligned.bam
    else
      ln -s aligned.bam ~{sample_id}.~{movie}.~{ref_name}.aligned.bam
      ln -s aligned.bam.bai ~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai
    fi

    wait ${BAM_STATS_PID}
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