version 1.0

import "../structs.wdl"

task sawfish_discover {
  meta {
    description: "Discover structural variant signatures with sawfish."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
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
    exclude_bed: {
      name: "Regions to exclude from CNV calls"
    }
    exclude_bed_index: {
      name: "Regions to exclude from CNV calls (index)"
    }
    expected_male_bed: {
      name: "Expected CN BED for sample with XY karyotype"
    }
    expected_female_bed: {
      name: "Expected CN BED for sample with XX karyotype"
    }
    small_variant_vcf: {
      name: "Small variant VCF"
    }
    small_variant_vcf_index: {
      name: "Small variant VCF index"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    discover_tar: {
      name: "Tarballed output of sawfish discover"
    }
    msg: {
      name: "Array of messages"
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

    File small_variant_vcf
    File small_variant_vcf_index

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  File expected_bed = if select_first([sex, "FEMALE"]) == "MALE" then expected_male_bed else expected_female_bed

  Int threads   = 16
  Int mem_gb    = threads * 8
  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    if [ "~{defined(sex)}" != "true" ]; then
      echo "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for sawfish."} >> messages.txt
    fi

    # sawfish stores relative filepaths of input files in the output directory
    # symlink input files into the working directory
    for i in ~{aligned_bam} ~{aligned_bam_index} ~{ref_fasta} ~{ref_index}; do
      ln --symbolic $i .
    done

    sawfish --version

    sawfish discover \
      --threads ~{threads} \
      --disable-path-canonicalization \
      --ref ~{basename(ref_fasta)} \
      --bam ~{basename(aligned_bam)} \
      --expected-cn ~{expected_bed} \
      --cnv-excluded-regions ~{exclude_bed} \
      --maf ~{small_variant_vcf} \
      --output-dir ~{out_prefix}

    tar --create --verbose --file ~{out_prefix}.tar ~{out_prefix}
    rm --recursive --force --verbose ~{out_prefix}
  >>>

  output {
    File discover_tar = "~{out_prefix}.tar"
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/sawfish@sha256:18ba096219fea38d6b32f5706fb794a05cc5d1d6cc16e2a09e3a13d62d8181d4"
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

task sawfish_call {
  meta {
    description: "Call structural variants from signatures with sawfish call."
  }

  parameter_meta {
    sample_ids: {
      name: "Sample IDs"
    }
    discover_tars: {
      name: "Tarballed output of sawfish discover"
    }
    aligned_bams: {
      name: "Aligned BAMs"
    }
    aligned_bam_indices: {
      name: "Aligned BAM indices"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    out_prefix: {
      name: "Output prefix"
    }
    report_supporting_reads: {
      name: "Report supporting reads"
    }
    copynum_bedgraph_names: {
      name: "Copy number bedgraph output filenames"
    }
    depth_bw_names: {
      name: "Depth bigWig output filenames"
    }
    gc_bias_corrected_depth_bw_names: {
      name: "GC bias corrected depth bigWig output filenames"
    }
    maf_bw_names: {
      name: "MAF bigWig output filenames"
    }
    copynum_summary_names: {
      name: "Copy number summary JSON output filenames"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    vcf: {
      name: "Structural variant VCF"
    }
    vcf_index: {
      name: "Structural variant VCF index"
    }
    supporting_reads: {
      name: "Supporting reads JSON"
    }
    copynum_bedgraph: {
      name: "Copy number bedgraph"
    }
    depth_bw: {
      name: "Depth bigWig"
    }
    gc_bias_corrected_depth_bw: {
      name: "GC bias corrected depth bigWig"
    }
    maf_bw: {
      name: "MAF bigWig"
    }
    copynum_summary: {
      name: "Copy number summary JSON"
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
    Array[String] maf_bw_names
    Array[String] copynum_summary_names

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = threads * 2
  Int disk_size = ceil(size(aligned_bams, "GB") + size(ref_fasta, "GB") + (size(discover_tars, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    # sawfish stores relative filepaths of input files in the output directory
    # symlink input files into the working directory
    for i in ~{sep=" " aligned_bams} ~{sep=" " aligned_bam_indices} ~{ref_fasta} ~{ref_index}; do
      ln --symbolic $i .
    done

    # unpack sawfish discover output
    # add the dir names to the SAMPLES list
    PREFIX="--sample "
    SAMPLES=()

    while read -r discover_tar || [[ -n "${discover_tar}" ]]; do
      sampledir=$(basename -s .tar "${discover_tar}")
      SAMPLES+=("$sampledir")
      tar --no-same-owner --extract --verbose --file "${discover_tar}"
    done < ~{write_lines(discover_tars)}

    sawfish --version

    # shellcheck disable=SC2068
    sawfish joint-call \
      --threads ~{threads} \
      ~{if report_supporting_reads then "--report-supporting-reads" else ""} \
      ${SAMPLES[@]/#/$PREFIX} \
      --output-dir ~{out_prefix}

    # sawshark annotation
    sawshark \
      --threads ~{(threads / 2)} \
      --vcf ~{out_prefix}/genotyped.sv.vcf.gz \
    | bcftools view - \
      --threads ~{(threads / 2) - 1} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz
    bcftools index --tbi \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{out_prefix}.vcf.gz

    # rename the output files to be more informative
    mv --verbose ~{out_prefix}/supporting_reads.json.gz ~{out_prefix}.supporting_reads.json.gz

    for sample_id in ~{sep=" " sample_ids}; do
      if [ "~{length(sample_ids)}" -gt 1 ]; then
        PREFIX="${sample_id}.~{out_prefix}"
      else
        PREFIX="~{out_prefix}"
      fi
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/copynum.bedgraph ${PREFIX}.copynum.bedgraph
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/depth.bw ${PREFIX}.depth.bw
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/gc_bias_corrected_depth.bw ${PREFIX}.gc_bias_corrected_depth.bw
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/maf.bw ${PREFIX}.maf.bw
      mv --verbose ~{out_prefix}/samples/sample????_${sample_id}/copynum.summary.json ${PREFIX}.copynum.summary.json
    done

    # shellcheck disable=SC2086,SC2048
    rm --recursive --force --verbose ${SAMPLES[*]}
  >>>

  output {
    File  vcf                              = "~{out_prefix}.vcf.gz"
    File  vcf_index                        = "~{out_prefix}.vcf.gz.tbi"
    File? supporting_reads                 = "~{out_prefix}.supporting_reads.json.gz"
    Array[File] copynum_bedgraph           = copynum_bedgraph_names
    Array[File] depth_bw                   = depth_bw_names
    Array[File] gc_bias_corrected_depth_bw = gc_bias_corrected_depth_bw_names
    Array[File] maf_bw                     = maf_bw_names
    Array[File] copynum_summary            = copynum_summary_names
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/sawfish@sha256:18ba096219fea38d6b32f5706fb794a05cc5d1d6cc16e2a09e3a13d62d8181d4"
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