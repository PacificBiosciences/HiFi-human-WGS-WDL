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
      haplotagged_bams: {
        description: "Haplotagged BAMs"
      },
      haplotagged_bam_indices: {
        description: "Indices for haplotagged BAMs"
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
    aligned_bams: {
      description: "BAMs"
    }
    aligned_bam_indices: {
      description: "BAM indices"
    }
    haplotagged_bam_names: {
      description: "Haplotagged BAM names"
    }
    haplotagged_bam_index_names: {
      description: "Haplotagged BAM index names"
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
    String sample_id
    Array[File] vcfs
    Array[File] vcf_indices
    Array[String] phased_vcf_names
    Array[String] phased_vcf_index_names
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    Array[String] haplotagged_bam_names
    Array[String] haplotagged_bam_index_names
    String ref_name
    File ref_fasta
    File ref_index
    String? preset
    Int? min_gq
    Int? min_mapq
    Boolean disable_global_realignment = false
    Boolean no_supplemental_joins = false
    Boolean phase_singletons = false
    Int threads = 16
    Int mem_gb = 96
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") + size(ref_fasta, "GB") + size(aligned_bams, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    BAM_PREFIX="--bam "
    BAMS=()
    for i in ~{sep=" " aligned_bams}; do
      ln --symbolic --verbose "${i}" .
      # shellcheck disable=SC2086
      BAMS+=("$(basename ${i})")
    done
    for i in ~{sep=" " aligned_bam_indices}; do
      ln --symbolic --verbose "${i}" .
    done

    VCF_PREFIX="--vcf "
    VCFS=()
    for i in ~{sep=" " vcfs}; do
      ln --symbolic --verbose "${i}" .
      # shellcheck disable=SC2086
      VCFS+=("$(basename ${i})")
    done
    for i in ~{sep=" " vcf_indices}; do
      ln --symbolic --verbose "${i}" .
    done

    hiphase --version

    # shellcheck disable=SC2068,SC2086
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
        then "--min-vcf-qual '" + min_gq + "'"
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
      ${VCFS[@]/#/$VCF_PREFIX} \
      ~{sep=" " prefix("--output-vcf ", phased_vcf_names)} \
      ${BAMS[@]/#/$BAM_PREFIX} \
      ~{sep=" " prefix("--output-bam ", haplotagged_bam_names)} \
      --reference "~{basename(ref_fasta)}" \
      --summary-file "~{sample_id}.~{ref_name}.hiphase.stats.tsv" \
      --blocks-file "~{sample_id}.~{ref_name}.hiphase.blocks.tsv" \
      --haplotag-file "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv"

    gzip "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv" &
    PID=$!

    # pull the phased basepairs and phase block N50
    cat << EOF > get_tsv_stats.py
    import sys, pandas as pd
    df = pd.read_csv('~{sample_id}.~{ref_name}.hiphase.stats.tsv', sep='\t')
    print(df[df['chromosome'] == 'all'][sys.argv[1]].values[0])
    EOF

    python3 get_tsv_stats.py basepairs_per_block_sum > phased_basepairs.txt
    python3 get_tsv_stats.py block_ng50 > phase_block_ng50.txt

    wait ${PID}
  >>>

  output {
    Array[File] phased_vcfs = phased_vcf_names
    Array[File] phased_vcf_indices = phased_vcf_index_names
    Array[File] haplotagged_bams = haplotagged_bam_names
    Array[File] haplotagged_bam_indices = haplotagged_bam_index_names
    File phase_stats = "~{sample_id}.~{ref_name}.hiphase.stats.tsv"
    File phase_blocks = "~{sample_id}.~{ref_name}.hiphase.blocks.tsv"
    File phase_haplotags = "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv.gz"
    String stat_phased_basepairs = read_string("phased_basepairs.txt")
    String stat_phase_block_ng50 = read_string("phase_block_ng50.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hiphase@sha256:41ebe22b55c66e2e78da2013f7fffaecc02a8b4e980400c3ea8d03c87330522e"  # 1.7.0_build2
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

