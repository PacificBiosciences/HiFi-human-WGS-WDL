version 1.0

import "../../structs.wdl"

workflow pbmm2 {
  meta {
    description: "Align HiFi reads to a reference genome with pbmm2"
    category: "Utility - Alignment"
    outputs: {
      aligned_bams: {
        description: "Array of aligned BAMs"
      },
      aligned_bam_indices: {
        description: "Array of aligned BAM indices"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    bam: {
      description: "Unaligned BAM"
    }
    use_alignment_chunking: {
      description: "Whether to split the input BAM into chunks for parallel alignment"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_name: {
      description: "Reference name"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    Boolean use_alignment_chunking = true
    File ref_fasta
    String ref_name
    RuntimeAttributes default_runtime_attributes
  }

  call create_pbmm2_index { input:
    ref_fasta = ref_fasta,
    runtime_attributes = default_runtime_attributes
  }

  call index_input_bam { input:
    bam = bam,
    use_alignment_chunking = use_alignment_chunking,
    runtime_attributes = default_runtime_attributes
  }

  scatter (chunk_index in if (index_input_bam.num_chunks > 0)
    then range(index_input_bam.num_chunks)
    else [
      -1
    ]
  ) {
    call pbmm2_align_wgs { input:
      sample_id = sample_id,
      bam = bam,
      pbindex = index_input_bam.pbindex,
      pbmm2_index = create_pbmm2_index.index,
      ref_name = ref_name,
      chunk_index = chunk_index + 1,  # wdl is 0-indexed, pbmm2 chunking is 1-indexed
      num_chunks = index_input_bam.num_chunks,
      runtime_attributes = default_runtime_attributes
    }
  }

  output {
    Array[File] aligned_bams = pbmm2_align_wgs.aligned_bam
    Array[File] aligned_bam_indices = pbmm2_align_wgs.aligned_bam_index
    Array[String] msg = index_input_bam.msg
  }
}

task create_pbmm2_index {
  meta {
    description: "Create a pbmm2 index for a reference genome"
    outputs: {
      index: {
        description: "pbmm2 reference index"
      }
    }
  }

  parameter_meta {
    ref_fasta: {
      description: "Reference FASTA"
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
    File ref_fasta
    Int threads = 8
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{ref_fasta}" .

    pbmm2 --version

    pbmm2 index \
      --num-threads ~{threads} \
      --log-level INFO \
      --preset HIFI \
      "~{basename(ref_fasta)}" \
      "~{basename(ref_fasta, ".fasta")}.mmi"
  >>>

  output {
    File index = "~{basename(ref_fasta, ".fasta")}.mmi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:86e39aa67fa5d385769d5f119739e8811ce550163e3b9dfc42bd58d1fecdf3a8"  # 26.2.99_build1
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

task index_input_bam {
  meta {
    description: "Split HiFi uBAM into chunks of a max size"
    outputs: {
      pbindex: {
        description: "pbindex for the input BAM"
      },
      num_chunks: {
        description: "Number of BAM chunks"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    bam: {
      description: "HiFi reads (BAM)"
    }
    use_alignment_chunking: {
      description: "Whether to split the input BAM into chunks for parallel alignment"
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
    File bam
    Boolean use_alignment_chunking
    Int threads = 16
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(bam, "GB") * 3 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{bam}" .

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
    print(json.dumps(check_bam_file('~{basename(bam)}', 10000)))
    EOF

    read -r kinetics base_modification aligned <<< "$(python3 ./detect_bam_tags.py | jq -r '. | [.kinetics, .base_modification, .aligned] | @tsv')"

    if [ "${aligned}" = true ]; then
      echo "Input ~{basename(bam)} is already aligned.  Alignments and haplotype tags will be stripped, and chunking will be disabled." >> messages.txt
    fi

    if [ "${base_modification}" = false ]; then
      echo "Input ~{basename(bam)} does not contain base modification tags.  5mCpG pileups will not be generated." >> messages.txt
    fi

    if [ "${kinetics}" = true ]; then
      echo "Input ~{basename(bam)} contains consensus kinetics tags. Kinetics will be stripped from the output." >> messages.txt
    fi

    echo 0 > chunks.txt  # sentinel value for no chunking

    # if chunking is desired, index the input BAM
    if [ "~{use_alignment_chunking}" = "true" ] && [ "${aligned}" = "false" ]; then
      pbindex --version
      pbindex --num-threads ~{threads} "~{basename(bam)}"

      echo 16 > chunks.txt
    fi
  >>>

  output {
    File? pbindex = "~{basename(bam)}.pbi"
    Int num_chunks = read_int("chunks.txt")
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbtk@sha256:f27bafa0ae6ffff6170d45a86a5677406320c7866760391f89ebe900a74d2039"  # 3.5.0_build2
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

task pbmm2_align_wgs {
  meta {
    description: "Align HiFi reads to a reference genome with pbmm2"
    outputs: {
      aligned_bam: {
        description: "Aligned BAM"
      },
      aligned_bam_index: {
        description: "Aligned BAM index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    bam: {
      description: "Unaligned BAM"
    }
    pbindex: {
      description: "pbindex for the input BAM"
    }
    pbmm2_index: {
      description: "pbmm2 reference index"
    }
    ref_name: {
      description: "Reference name"
    }
    strip_kinetics: {
      description: "Strip kinetics tags"
    }
    keep_unmapped: {
      description: "Keep unmapped reads"
    }
    min_length: {
      description: "Minimum mapped read length in basepairs"
    }
    threads: {
      description: "Number of threads to use"
    }
    mem_gb: {
      description: "Memory to use in GiB"
    }
    chunk_index: {
      description: "Chunk index for parallel alignment"
    }
    num_chunks: {
      description: "Total number of chunks for parallel alignment"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File bam
    File? pbindex
    File pbmm2_index
    String ref_name
    Boolean strip_kinetics = true
    Boolean keep_unmapped = true
    Int min_length = 50
    Int threads = 32
    Int mem_gb = 64
    Int chunk_index = 0
    Int num_chunks = 1
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(bam, "GB") * 2 + size(pbmm2_index, "GB") + 70)

  Boolean use_chunking = chunk_index > 0
  String movie = basename(bam, ".bam")

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bam}" "~{pbmm2_index}" .
    ~{if (defined(pbindex))
      then "ln --symbolic --verbose ~{pbindex} ~{basename(bam)}.pbi"
      else ""
    }

    pbmm2 --version
    pbmm2 align \
      --num-threads ~{threads} \
      --sort-memory 4G \
      --preset HIFI \
      --sample "~{sample_id}" \
      --log-level INFO \
      --sort \
      --strip-tags HP,PS,PC \
      ~{true="--strip" false="" strip_kinetics} \
      ~{true="--unmapped" false="" keep_unmapped} \
      --min-length ~{min_length} \
      ~{if (use_chunking)
        then "--chunk '" + chunk_index + "/" + num_chunks + "' --chunk-mode scatter"
        else ""
      } \
      "~{basename(pbmm2_index)}" \
      "~{basename(bam)}" \
      "~{sample_id}.~{movie}.chunk_~{chunk_index}.~{ref_name}.aligned.bam"
  >>>

  output {
    File aligned_bam = "~{sample_id}.~{movie}.chunk_~{chunk_index}.~{ref_name}.aligned.bam"
    File aligned_bam_index = "~{sample_id}.~{movie}.chunk_~{chunk_index}.~{ref_name}.aligned.bam.bai"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:86e39aa67fa5d385769d5f119739e8811ce550163e3b9dfc42bd58d1fecdf3a8"  # 26.2.99_build1
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

