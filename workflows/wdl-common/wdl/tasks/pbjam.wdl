version 1.0

import "../structs.wdl"

task pbjam_bam_stats {
  meta {
    description: "Get read, mapping, and alignment statistics from a BAM file."
    outputs: {
      read_length_plot: {
        description: "Distribution of read lengths"
      },
      read_quality_plot: {
        description: "Distribution of read qualities"
      },
      mapq_distribution_plot: {
        description: "Distribution of mapping quality per alignment"
      },
      mg_distribution_plot: {
        description: "Distribution of gap-compressed identity per alignment"
      },
      stat_read_count: {
        description: "Number of reads"
      },
      stat_read_length_mean: {
        description: "Mean read length"
      },
      stat_read_length_median: {
        description: "Median read length"
      },
      stat_read_length_n50: {
        description: "Read length N50"
      },
      stat_read_quality_mean: {
        description: "Mean read quality"
      },
      stat_read_quality_median: {
        description: "Median read quality"
      },
      stat_mapped_read_count: {
        description: "Number of reads mapped to reference"
      },
      stat_mapped_read_percent: {
        description: "Percent of reads mapped to reference"
      },
      stat_gap_compressed_identity_mean: {
        description: "Mean gap-compressed identity"
      },
      stat_gap_compressed_identity_median: {
        description: "Median gap-compressed identity"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    ref_name: {
      description: "Reference name"
    }
    bam: {
      description: "Input BAM"
    }
    bam_index: {
      description: "Input BAM index"
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
    String ref_name
    File bam
    File bam_index
    Int threads = 8
    Int mem_gb = 8
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(bam, "GB") + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bam}" "~{bam_index}" .

    pbjam bam-stats \
      --threads ~{threads} \
      --include-unmapped \
      --input-bam "~{basename(bam)}" \
      --output-json "~{sample_id}.~{ref_name}.pbjam.json" \
      --plot-label "~{sample_id}.~{ref_name}" \
      --image-prefix "~{sample_id}.~{ref_name}"

    mv --verbose "~{sample_id}.~{ref_name}.read_length_distribution.png" "~{sample_id}.read_length_distribution.png"
    mv --verbose "~{sample_id}.~{ref_name}.read_quality_distribution.png" "~{sample_id}.read_quality_distribution.png"

    cat << EOF > get_value.py
    import json, decimal, sys
    data = json.load(open("~{sample_id}.~{ref_name}.pbjam.json"), parse_float=decimal.Decimal)
    print(data['combined_stats']['summary_stats'][sys.argv[1]])
    EOF

    python3 ./get_value.py read_count > read_count.txt
    python3 ./get_value.py read_length_mean > read_length_mean.txt
    python3 ./get_value.py read_length_median > read_length_median.txt
    python3 ./get_value.py read_length_n50 > read_length_n50.txt
    python3 ./get_value.py read_quality_mean > read_quality_mean.txt
    python3 ./get_value.py read_quality_median > read_quality_median.txt
    python3 ./get_value.py mapped_read_count > mapped_read_count.txt
    python3 ./get_value.py mapped_read_percentage > mapped_read_percent.txt
    python3 ./get_value.py gap_compress_identity_mean > mg_mean.txt
    python3 ./get_value.py gap_compress_identity_median > mg_median.txt
  >>>

  output {
    File read_length_plot = "~{sample_id}.read_length_distribution.png"
    File read_quality_plot = "~{sample_id}.read_quality_distribution.png"
    File mapq_distribution_plot = "~{sample_id}.~{ref_name}.mapq_distribution.png"
    File mg_distribution_plot = "~{sample_id}.~{ref_name}.mg_distribution.png"
    String stat_read_count = read_string("read_count.txt")
    String stat_read_length_mean = read_string("read_length_mean.txt")
    String stat_read_length_median = read_string("read_length_median.txt")
    String stat_read_length_n50 = read_string("read_length_n50.txt")
    String stat_read_quality_mean = read_string("read_quality_mean.txt")
    String stat_read_quality_median = read_string("read_quality_median.txt")
    String stat_mapped_read_count = read_string("mapped_read_count.txt")
    String stat_mapped_read_percent = read_string("mapped_read_percent.txt")
    String stat_gap_compressed_identity_mean = read_string("mg_mean.txt")
    String stat_gap_compressed_identity_median = read_string("mg_median.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbjam@sha256:1b90537ed6683ac7b89002c0f59e940b4ea46c600d292904f549af1ae2760bba"  # 0.4.1_build2
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

