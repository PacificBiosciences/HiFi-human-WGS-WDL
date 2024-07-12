version 1.0

import "../structs.wdl"

task mosdepth {
  meta {
    description: "Calculate coverage statistics using mosdepth"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    ref_name: {
      name: "Reference name"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    infer_sex: {
      name: "Infer the sex of human samples based on chrY depth"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    summary: {
      name: "Summary depth statistics"
    }
    region_bed: {
      name: "Depth BED"
    }
    inferred_sex: {
      name: "Sex inferred from chrY depth"
    }
    stat_mean_depth: {
      name: "Mean depth"
    }
  }

  input {
    String sample_id
    String ref_name

    File aligned_bam
    File aligned_bam_index

    Boolean infer_sex = false

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 4
  Int disk_size = ceil(size(aligned_bam, "GB") + 20)

  Float max_norm_female_chrY_depth = 0.1

  String out_prefix = basename(aligned_bam, ".bam")

  command <<<
    set -euo pipefail

    mosdepth --version

    mosdepth \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --by 500 \
      --no-per-base \
      --use-median \
      ~{out_prefix} \
      ~{aligned_bam}

    # normalize output names
    if [ ! -f ~{sample_id}.~{ref_name}.mosdepth.summary.txt ]; then
      mv ~{out_prefix}.mosdepth.summary.txt ~{sample_id}.~{ref_name}.mosdepth.summary.txt
      mv ~{out_prefix}.regions.bed.gz ~{sample_id}.~{ref_name}.regions.bed.gz
    fi

    cat << EOF > get_mean_depth.py
    import pandas as pd
    df = pd.read_csv('~{sample_id}.~{ref_name}.mosdepth.summary.txt', sep='\t')
    print(df[df['chrom'] == 'total']['mean'].values[0])
    EOF

    python3 ./get_mean_depth.py > mean_depth.txt || echo "0.0" > mean_depth.txt

    cat << EOF > infer_chrY.py
    import pandas as pd
    df = pd.read_csv('~{sample_id}.~{ref_name}.mosdepth.summary.txt', sep='\t')
    chrA_depth = df[df['chrom'].str.match(r'^(chr)?\d{1,2}$')]['mean'].mean()
    chrY_depth = df[df['chrom'].str.match(r'^(chr)?Y$')]['mean'].mean()
    if chrA_depth == 0:
      print("")
    else:
      print("MALE" if chrY_depth/chrA_depth > float(~{max_norm_female_chrY_depth}) else "FEMALE")
    EOF

    python3 ./infer_chrY.py > inferred_sex.txt || echo "" > inferred_sex.txt
  >>>

  output {
    File   summary         = "~{sample_id}.~{ref_name}.mosdepth.summary.txt"
    File   region_bed      = "~{sample_id}.~{ref_name}.regions.bed.gz"
    String inferred_sex    = if (infer_sex) then read_string("inferred_sex.txt") else ""
    String stat_mean_depth = read_string("mean_depth.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:f715c11100e9bb3562cce1c5e23a185cfcc92a6fec412b16c30c0250496cc0d1"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " LOCAL"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
