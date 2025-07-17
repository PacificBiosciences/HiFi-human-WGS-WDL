version 1.0

import "../../structs.wdl"

workflow pbmm2 {
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
    max_reads_per_chunk: {
      name: "Maximum reads per alignment chunk"
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
    default_runtime_attributes: {
      name: "Default runtime attribute structure"
    }
    aligned_bams: {
      name: "Array of aligned BAMs"
    }
    aligned_bam_indices: {
      name: "Array of aligned BAM indices"
    }
    msg: {
      name: "Array of messages"
    }
  }

  input {
    String sample_id
    File bam

    Int max_reads_per_chunk

    File ref_fasta
    File ref_index
    String ref_name

    RuntimeAttributes default_runtime_attributes
  }

  call split_input_bam {
    input:
      bam = bam,
      max_reads_per_chunk = max_reads_per_chunk,
      runtime_attributes = default_runtime_attributes
  }

  scatter (chunk in if (length(split_input_bam.chunks) > 0) then split_input_bam.chunks else [bam] ) {
    call pbmm2_align_wgs {
      input:
        sample_id = sample_id,
        bam = chunk,
        ref_fasta = ref_fasta,
        ref_index = ref_index,
        ref_name = ref_name,
        runtime_attributes = default_runtime_attributes
    }
  }

  output {
    Array[File] aligned_bams        = pbmm2_align_wgs.aligned_bam
    Array[File] aligned_bam_indices = pbmm2_align_wgs.aligned_bam_index
    Array[String] msg               = split_input_bam.msg
  }
}

task split_input_bam {
  meta {
    description: "Split HiFi uBAM into chunks of a max size."
  }

  parameter_meta {
    bam: {
      name: "HiFi reads (BAM)"
    }
    max_reads_per_chunk: {
      name: "Maximum reads per chunk"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    chunks: {
      name: "Split BAMs"
    }
    msg: {
      name: "Array of messages"
    }
  }

  input {
    File bam
    Int max_reads_per_chunk

    RuntimeAttributes runtime_attributes
  }

  String movie = basename(bam, ".bam")

  Int threads   = 16
  Int mem_gb    = 32
  Int disk_size = ceil(size(bam, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    INBAM="~{bam}"

    # Check for presence of alignment, basemod, and kinetics tags
    cat << EOF > detect_bam_tags.py
    import json, pysam
    def check_bam_file(bam_file_path, n_records):
      output = dict()
      save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
      with pysam.AlignmentFile(bam_file_path, 'rb', check_sq=False) as bam_file:
        pysam.set_verbosity(save)  # restore warnings
        aligned = bool(bam_file.nreferences)
        unique_tags = set()
        for i, record in enumerate(bam_file):
          if i >= n_records: break
          unique_tags.update(tag[0] for tag in record.tags)
      output['kinetics'] = bool(unique_tags & {'fi', 'ri', 'fp', 'rp', 'ip', 'pw'})
      output['base_modification'] = bool(unique_tags & {'MM', 'ML', 'Mm', 'Ml'})
      output['aligned'] = aligned
      return output
    print(json.dumps(check_bam_file('~{bam}', 10000)))
    EOF

    read -r kinetics base_modification aligned <<< "$(python3 ./detect_bam_tags.py | jq -r '. | [.kinetics, .base_modification, .aligned] | @tsv')"

    if [ "$aligned" = true ]; then
      echo "Input ~{basename(bam)} is already aligned.  Alignments and haplotype tags will be stripped." >> messages.txt
    fi

    if [ "$base_modification" = false ]; then
      echo "Input ~{basename(bam)} does not contain base modification tags.  5mCpG pileups will not be generated." >> messages.txt
    fi

    if [ "$kinetics" = true ]; then
      echo "Input ~{basename(bam)} contains consensus kinetics tags. Kinetics will be stripped from the output." >> messages.txt
    fi

    # reset BAM and strip kinetics/haplotype tags if present
    if [ "$aligned" = true ] || [ "$kinetics" = true ]; then
      samtools --version
      samtools reset \
        ~{if threads > 1 then "-@ " + (threads - 1) else ""} \
        --remove-tag fi,ri,fp,rp,ip,pw,HP,PS,PC \
        -o ~{movie}.reset.bam \
        ~{bam}
        INBAM=~{movie}.reset.bam
    fi

    # if chunking is desired, index the input BAM and list ZMWs
    if [ "~{max_reads_per_chunk}" -gt "1" ]; then
      pbindex --version
      pbindex --num-threads ~{threads} $INBAM

      zmwfilter --version
      zmwfilter --num-threads ~{threads} --show-all $INBAM > ~{movie}.zmws.txt

      read -r n_records <<< "$(wc -l ~{movie}.zmws.txt | cut -f1 -d' ')"

      # if the number of ZMWs is greater than the chunk size, split the input BAM
      if [ "$n_records" -gt "~{max_reads_per_chunk}" ]; then
        split --version
        split \
          --verbose \
          --lines=~{max_reads_per_chunk} \
          --numeric-suffixes \
          ~{movie}.zmws.txt \
          chunk_

        parallel --version
        # shellcheck disable=SC2012
        ls chunk_* | parallel -v -j ~{threads} \
          zmwfilter --num-threads 1 --include {} $INBAM ~{movie}.{}.bam

        # if the input BAM was reset, remove so that it is not included in the output
        rm --force --verbose ~{movie}.reset.bam
      fi
    fi
  >>>

  output {
    Array[File] chunks = glob("~{movie}.*.bam")
    Array[String] msg  = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbtk@sha256:67cd438ed9f343f90f058108170ddbff8fb1d9b5c193f4016be42b737ee2e73c"
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

  Int threads   = 32
  Int mem_gb    = ceil(threads * 4)
  Int disk_size = ceil(size(bam, "GB") * 2 + size(ref_fasta, "GB") + 70)

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
      ~{true='--strip' false='' strip_kinetics} \
      --unmapped \
      ~{ref_fasta} \
      ~{bam} \
      aligned.bam

    mv --verbose aligned.bam ~{sample_id}.~{movie}.~{ref_name}.aligned.bam
    mv --verbose aligned.bam.bai ~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai
  >>>

  output {
    File aligned_bam       = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
    File aligned_bam_index = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:5f3f4d1f5dbea5cd4c388ee26b2fecbbb7dbcef449343633e039dca3d3725859"
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