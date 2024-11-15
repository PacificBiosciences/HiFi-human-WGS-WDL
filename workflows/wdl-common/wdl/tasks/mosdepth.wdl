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
    region_bed_index: {
      name: "Depth BED index"
    }
    depth_distribution_plot: {
      name: "Depth distribution plot"
    }
    stat_mean_depth: {
      name: "Mean depth"
    }
    inferred_sex: {
      name: "Sex inferred from chrY depth"
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
      mv --verbose ~{out_prefix}.mosdepth.summary.txt ~{sample_id}.~{ref_name}.mosdepth.summary.txt
      mv --verbose ~{out_prefix}.regions.bed.gz ~{sample_id}.~{ref_name}.regions.bed.gz
      mv --verbose ~{out_prefix}.regions.bed.gz.csi ~{sample_id}.~{ref_name}.regions.bed.gz.csi
    fi

    # plot depth distribution
    cat << EOF > plot_depth_distribution.py
    import pandas as pd, seaborn as sns, matplotlib.pyplot as plt, numpy as np
    sns.set_theme(style='darkgrid')
    df = pd.read_csv(
      '~{sample_id}.~{ref_name}.regions.bed.gz',
      sep='\t',
      names=['chr', 'start', 'end', 'depth'],
      compression='gzip',
    )
    xmax = int(2*df[df['depth'] > 0]['depth'].mode().values[0])  # 2x non-zero mode
    df = df[(df['depth'] <= xmax)]
    fig, axs = plt.subplots(2, 1, figsize=(8,6))
    sns.histplot(df, x='depth', bins=xmax-1, stat='proportion', ax=axs[0])
    axs[0].set_xlim(0, xmax)
    axs[0].set_xlabel('')
    axs[0].set_xticklabels([])
    sns.ecdfplot(df, x='depth', complementary=True, ax=axs[1])
    axs[1].set_xlim(0, xmax)
    axs[1].set_yticks(np.arange(0,1.1,0.1))
    axs[1].set_ylabel('Proportion ≥ depth')
    plt.suptitle('~{sample_id}.~{ref_name}\nAligned depth distribution')
    fig.tight_layout()
    plt.savefig('~{sample_id}.~{ref_name}.depth_distribution.png')
    EOF

    python3 ./plot_depth_distribution.py

    # get the mean depth
    cat << EOF > get_mean_depth.py
    import pandas as pd
    df = pd.read_csv('~{sample_id}.~{ref_name}.mosdepth.summary.txt', sep='\t')
    print(df[df['chrom'] == 'total']['mean'].values[0])
    EOF

    python3 ./get_mean_depth.py > mean_depth.txt || echo "0.0" > mean_depth.txt

    # infer the sex from normalized chrY depth
    cat << EOF > infer_chrY.py
    import pandas as pd
    df = pd.read_csv('~{sample_id}.~{ref_name}.mosdepth.summary.txt', sep='\t')
    chrA_depth = df[df['chrom'].str.match(r'^(chr)?\d{1,2}$')]['mean'].mean()
    chrY_depth = df[df['chrom'].str.match(r'^(chr)?Y$')]['mean'].mean()
    if chrA_depth == 0:
      print('')
    else:
      print('MALE' if chrY_depth/chrA_depth > float(~{max_norm_female_chrY_depth}) else 'FEMALE')
    EOF

    python3 ./infer_chrY.py > inferred_sex.txt || echo "" > inferred_sex.txt
  >>>

  output {
    File   summary                 = "~{sample_id}.~{ref_name}.mosdepth.summary.txt"
    File   region_bed              = "~{sample_id}.~{ref_name}.regions.bed.gz"
    File   region_bed_index        = "~{sample_id}.~{ref_name}.regions.bed.gz.csi"
    File   depth_distribution_plot = "~{sample_id}.~{ref_name}.depth_distribution.png"
    String stat_mean_depth         = read_string("mean_depth.txt")
    String inferred_sex            = if (infer_sex) then read_string("inferred_sex.txt") else ""
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:63f7a5d1a4a17b71e66d755d3301a951e50f6b63777d34dab3ee9e182fd7acb1"
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
