version 1.0

import "../structs.wdl"

task sawfish_discover {
  meta {
    description: "Discover structural variant signatures with Sawfish"
    outputs: {
      discover_tar: {
        description: "Tarballed intermediate output of sawfish discover"
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
    sex: {
      description: "Sample sex",
      choices: [
        "MALE",
        "FEMALE"
      ]
    }
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
    exclude_bed: {
      description: "Regions to exclude from CNV calls"
    }
    exclude_bed_index: {
      description: "Regions to exclude from CNV calls (index)"
    }
    expected_male_bed: {
      description: "Expected CN BED for sample with XY karyotype"
    }
    expected_female_bed: {
      description: "Expected CN BED for sample with XX karyotype"
    }
    out_prefix: {
      description: "Output prefix"
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
    String? sex
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    File exclude_bed
    File exclude_bed_index
    File expected_male_bed
    File expected_female_bed
    String out_prefix
    Int threads = 16
    Int mem_gb = 128
    RuntimeAttributes runtime_attributes
  }

  File expected_bed = if select_first([
    sex,
    "FEMALE"
  ]) == "MALE"
    then expected_male_bed
    else expected_female_bed

  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    if [ "~{defined(sex)}" != "true" ]; then
      echo "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for sawfish." >> messages.txt
    fi

    ln --symbolic --verbose "~{aligned_bam}" "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .
    ln --symbolic --verbose "~{exclude_bed}" "~{exclude_bed_index}" .
    ln --symbolic --verbose "~{expected_bed}" .

    sawfish --version
    sawfish discover \
      --threads ~{threads} \
      --disable-path-canonicalization \
      --ref "~{basename(ref_fasta)}" \
      --bam "~{basename(aligned_bam)}" \
      --expected-cn "~{basename(expected_bed)}" \
      --cnv-excluded-regions "~{basename(exclude_bed)}" \
      --output-dir "~{out_prefix}"

    tar --create --verbose --file "~{out_prefix}.tar" "~{out_prefix}"
    rm --recursive --force --verbose "~{out_prefix}"
  >>>

  output {
    File discover_tar = "~{out_prefix}.tar"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/sawfish@sha256:5cf8f02790ecb89e652885c57f103fc9597c41f75ee71be2c05510cbbdf68f59"  # 2.2.1_build2
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

task sawfish_call {
  meta {
    description: "Call structural variants from signatures with sawfish call"
    outputs: {
      vcf: {
        description: "Structural variant VCF"
      },
      vcf_index: {
        description: "Structural variant VCF index"
      },
      supporting_reads: {
        description: "Supporting reads for structural variants"
      },
      copynum_bedgraph: {
        description: "CNV copy number BEDGraph"
      },
      depth_bw: {
        description: "CNV depth BigWig"
      },
      gc_bias_corrected_depth_bw: {
        description: "CNV GC-bias corrected depth BigWig"
      },
      copynum_summary: {
        description: "CNV copy number summary JSON"
      }
    }
  }

  parameter_meta {
    sample_ids: {
      description: "Sample IDs"
    }
    discover_tars: {
      description: "Tarballed intermediate output of sawfish discover"
    }
    aligned_bams: {
      description: "Aligned BAMs"
    }
    aligned_bam_indices: {
      description: "Aligned BAM indices"
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
    report_supporting_reads: {
      description: "Report supporting reads"
    }
    copynum_bedgraph_names: {
      description: "CNV copy number BEDGraph output filenames"
    }
    depth_bw_names: {
      description: "CNV depth BigWig output filenames"
    }
    gc_bias_corrected_depth_bw_names: {
      description: "CNV GC-bias corrected depth BigWig output filenames"
    }
    copynum_summary_names: {
      description: "CNV copy number summary JSON output filenames"
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
    Array[String] sample_ids
    Array[File] discover_tars
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    File ref_fasta
    File ref_index
    String out_prefix
    Boolean report_supporting_reads = true
    Array[String] copynum_bedgraph_names
    Array[String] depth_bw_names
    Array[String] gc_bias_corrected_depth_bw_names
    Array[String] copynum_summary_names
    Int threads = 16
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(aligned_bams, "GB") + size(ref_fasta, "GB") + (size(discover_tars, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    # sawfish stores relative filepaths of input files in the output directory
    # symlink input files into the working directory
    for i in ~{sep=" " aligned_bams} ~{sep=" " aligned_bam_indices}; do
      ln --symbolic --verbose "${i}" .
    done
    ln --symbolic --verbose "~{ref_fasta}" "~{ref_index}" .

    # unpack sawfish discover output
    # add the dir names to the SAMPLES list
    PREFIX="--sample "
    SAMPLES=()

    while read -r discover_tar || [[ -n "${discover_tar}" ]]; do
      sampledir=$(basename -s .tar "${discover_tar}")
      SAMPLES+=("${sampledir}")
      tar --no-same-owner --extract --verbose --file "${discover_tar}"
    done < "~{write_lines(discover_tars)}"

    sawfish --version

    # shellcheck disable=SC2068
    sawfish joint-call \
      --threads ~{threads} \
      ~{if report_supporting_reads
        then "--report-supporting-reads"
        else ""
      } \
      ${SAMPLES[@]/#/$PREFIX} \
      --output-dir "~{out_prefix}"

    # sawshark annotation
    sawshark \
      ~{if threads > 1
        then "--threads '" + ceil(threads / 2) + "'"
        else ""
      } \
      --vcf "~{out_prefix}/genotyped.sv.vcf.gz" \
    | bcftools view - \
      ~{if (threads - ceil(threads / 2)) > 1
        then "--threads '" + (threads - ceil(threads / 2) - 1) + "'"
        else ""
      } \
      --output-type z \
      --write-index=tbi \
      --output "~{out_prefix}.vcf.gz"

    # rename the output files to be more informative
    # sawfish only writes supporting_reads.json.gz when
    # --report-supporting-reads was passed above.
    ~{if report_supporting_reads
      then "mv --verbose '" + out_prefix + "/supporting_reads.json.gz' '" + out_prefix + ".supporting_reads.json.gz'"
      else ""
    }

    for sample_id in ~{sep=" " sample_ids}; do
      if [ "~{length(sample_ids)}" -gt 1 ]; then
        PREFIX="${sample_id}.~{out_prefix}"
      else
        PREFIX="~{out_prefix}"
      fi
      # shellcheck disable=SC2086
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/copynum.bedgraph "${PREFIX}.copynum.bedgraph"
      # shellcheck disable=SC2086
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/depth.bw "${PREFIX}.depth.bw"
      # shellcheck disable=SC2086
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/gc_bias_corrected_depth.bw "${PREFIX}.gc_bias_corrected_depth.bw"
      # shellcheck disable=SC2086
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/copynum.summary.json "${PREFIX}.copynum.summary.json"
    done

    # shellcheck disable=SC2086,SC2048
    rm --recursive --force --verbose ${SAMPLES[*]}
  >>>

  output {
    File vcf = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File? supporting_reads = "~{out_prefix}.supporting_reads.json.gz"
    Array[File] copynum_bedgraph = copynum_bedgraph_names
    Array[File] depth_bw = depth_bw_names
    Array[File] gc_bias_corrected_depth_bw = gc_bias_corrected_depth_bw_names
    Array[File] copynum_summary = copynum_summary_names
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/sawfish@sha256:5cf8f02790ecb89e652885c57f103fc9597c41f75ee71be2c05510cbbdf68f59"  # 2.2.1_build2
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

