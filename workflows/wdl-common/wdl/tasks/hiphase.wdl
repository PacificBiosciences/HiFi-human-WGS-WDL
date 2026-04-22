version 1.0

import "../structs.wdl"

task hiphase {
  meta {
    description: "Phases VCFs and haplotags BAMs with HiPhase"
    outputs: {
      phased_vcfs: {
        description: "Phased VCFs"
      },
      phased_vcf_indices: {
        description: "Indices for phased VCFs"
      },
      haplotagged_bam: {
        description: "Haplotagged BAM"
      },
      haplotagged_bam_index: {
        description: "Index for haplotagged BAM"
      },
      phase_stats: {
        description: "Phasing statistics"
      },
      phase_blocks: {
        description: "Phase blocks"
      },
      phase_haplotags: {
        description: "Per-read phase assignment"
      },
      stat_phased_basepairs: {
        description: "Number of basepairs within phase blocks"
      },
      stat_phase_block_ng50: {
        description: "Phase block NG50"
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
    vcfs: {
      description: "VCFs"
    }
    vcf_indices: {
      description: "VCF indices"
    }
    phased_vcf_names: {
      description: "Phased VCF names"
    }
    phased_vcf_index_names: {
      description: "Phased VCF index names"
    }
    aligned_bam: {
      description: "BAM"
    }
    aligned_bam_index: {
      description: "BAM index"
    }
    ref_name: {
      description: "Reference name"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    preset: {
      description: "Preset configuration to use",
      help: "If unset, default parameters will be used"
    }
    min_gq: {
      description: "Minimum genotype quality to consider a variant for phasing"
    }
    min_mapq: {
      description: "Minimum mapping quality to consider a read for phasing"
    }
    disable_global_realignment: {
      description: "If true, disable global realignment of reads"
    }
    no_supplemental_joins: {
      description: "If true, do not use supplemental alignments to join phase blocks"
    }
    phase_singletons: {
      description: "If true, attempt to phase variants that are not connected to any other variants"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
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
    String? preset
    Int? min_gq
    Int? min_mapq
    Boolean disable_global_realignment = false
    Boolean no_supplemental_joins = false
    Boolean phase_singletons = false
    RuntimeAttributes runtime_attributes
  }

  Int threads = 16
  Int mem_gb = threads * 6
  Int disk_size = ceil(size(vcfs, "GB") + size(ref_fasta, "GB") + size(aligned_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{aligned_bam}" .
    ln --symbolic --verbose "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    hiphase --version

    # shellcheck disable=SC2086
    hiphase --threads ~{threads} \
      ~{if defined(preset)
        then "--preset '" + preset + "'"
        else ""
      } \
      ~{if defined(min_mapq)
        then "--min-mapq '" + min_mapq + "'"
        else ""
      } \
      ~{if defined(min_gq)
        then "--min-gq '" + min_gq + "'"
        else ""
      } \
      ~{if no_supplemental_joins
        then "--no-supplemental-joins"
        else ""
      } \
      ~{if phase_singletons
        then "--phase-singletons"
        else ""
      } \
      ~{if disable_global_realignment
        then "--disable-global-realignment"
        else ""
      } \
      --sample-name "~{sample_id}" \
      ~{sep=" " prefix("--vcf ", vcfs)} \
      ~{sep=" " prefix("--output-vcf ", phased_vcf_names)} \
      --bam "~{basename(aligned_bam)}" \
      --output-bam "~{sample_id}.~{ref_name}.haplotagged.bam" \
      --reference "~{basename(ref_fasta)}" \
      --summary-file "~{sample_id}.~{ref_name}.hiphase.stats.tsv" \
      --blocks-file "~{sample_id}.~{ref_name}.hiphase.blocks.tsv" \
      --haplotag-file "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv"

    gzip "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv"

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
    Array[File] phased_vcfs = phased_vcf_names
    Array[File] phased_vcf_indices = phased_vcf_index_names
    File haplotagged_bam = "~{sample_id}.~{ref_name}.haplotagged.bam"
    File haplotagged_bam_index = "~{sample_id}.~{ref_name}.haplotagged.bam.bai"
    File phase_stats = "~{sample_id}.~{ref_name}.hiphase.stats.tsv"
    File phase_blocks = "~{sample_id}.~{ref_name}.hiphase.blocks.tsv"
    File phase_haplotags = "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv.gz"
    String stat_phased_basepairs = read_string("phased_basepairs.txt")
    String stat_phase_block_ng50 = read_string("phase_block_ng50.txt")
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hiphase@sha256:2c54932e4992d5fcee03395f271a54daa69d27294d11cfa5a7d05964fbfb1ca6"  # 1.6.0_build1
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
