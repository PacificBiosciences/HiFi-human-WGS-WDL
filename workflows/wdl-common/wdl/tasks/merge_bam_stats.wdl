version 1.0

import "../structs.wdl"

task merge_bam_stats {
  meta {
    description: "Merge BAM stats files and create read length and quality histograms"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    bam_stats: {
      name: "BAM Stats"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    read_length_histogram: {
      name: "Read Length Summary"
    }
    read_quality_histogram: {
      name: "Read Quality Summary"
    }
    read_length_and_quality: {
      name: "Read Length and Quality"
    }
    stat_num_reads: {
      name: "Number of Reads"
    }
    stat_read_length_mean: {
      name: "Mean Read Length"
    }
    stat_read_length_median: {
      name: "Median Read Length"
    }
    stat_read_quality_mean: {
      name: "Mean Quality"
    }
    stat_read_quality_median: {
      name: "Median Quality"
    }
  }

  input {
    String sample_id
    Array[File] bam_stats

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = 10

  command <<<
    zcat ~{sep=" " bam_stats} > ~{sample_id}.read_length_and_quality.tsv

    cat << EOF > bin_length.py
    import sys, pandas as pd
    df = pd.read_csv(sys.stdin, sep='\t', header=None)
    # bin by read length, 0-40000, by 1000
    bins = pd.interval_range(start=0, end=40000, freq=1000, closed='left')
    df[0] = pd.cut(df[2], bins=bins, labels=bins.left)
    grouped = df.groupby(0)[2].agg(['count', 'sum'])
    # print the count and sum of read lengths for each bin
    for index, row in grouped.iterrows():
      print(f'{index.left}\t{row["count"]}\t{row["sum"]}')
    EOF

    python3 ./bin_length.py \
      < ~{sample_id}.read_length_and_quality.tsv \
      > ~{sample_id}.read_length_histogram.tsv

    cat << EOF > plot_length.py
    import sys, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
    sns.set_theme(style='darkgrid')
    df = pd.read_csv(sys.stdin, sep='\t', header=None, names=['read_length', 'count', 'base_pairs'])
    fig, axs = plt.subplots(1, 1, figsize=(7, 5))
    sns.histplot(data=df, x='read_length', weights='count', bins=40, ax=axs)
    axs.set_xlabel('read length (bp)'); axs.set_ylabel('read count')
    axs.set_title('~{sample_id} read length histogram')
    plt.tight_layout()
    plt.savefig('~{sample_id}.read_length_histogram.png')
    EOF

    python3 ./plot_length.py < ~{sample_id}.read_length_histogram.tsv

    cat << EOF > bin_quality.py
    import sys, pandas as pd
    df = pd.read_csv(sys.stdin, sep='\t', header=None)
    # bin by quality score, 0-60
    bins = pd.interval_range(start=0, end=61, freq=1, closed='left')
    df[0] = pd.cut(df[3], bins=bins, labels=bins.left, include_lowest=True)
    grouped = df.groupby(0)[2].agg(['count', 'sum'])
    # print the count and sum of read lengths for each quality score
    for index, row in grouped.iterrows():
      print(f'{index.left}\t{row["count"]}\t{row["sum"]}')
    EOF

    python3 ./bin_quality.py \
      < ~{sample_id}.read_length_and_quality.tsv \
      > ~{sample_id}.read_quality_histogram.tsv

    cat << EOF > plot_quality.py
    import sys, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
    sns.set_theme(style='darkgrid')
    df = pd.read_csv(sys.stdin, sep='\t', header=None, names=['read_quality', 'count', 'base_pairs'])
    fig, axs = plt.subplots(1, 1, figsize=(7, 5))
    sns.histplot(data=df, x='read_quality', weights='count', bins=60, ax=axs)
    axs.set_xlabel('read quality'); axs.set_ylabel('read count')
    axs.set_title('~{sample_id} read quality histogram')
    plt.tight_layout()
    plt.savefig('~{sample_id}.read_quality_histogram.png')
    EOF

    python3 ./plot_quality.py < ~{sample_id}.read_quality_histogram.tsv

    cat << EOF > summary_stats.py
    import sys, pandas as pd
    df = pd.read_csv(sys.stdin, sep='\t', header=None)
    print(f'{len(df)}\t{df[2].mean().round(2)}\t{df[2].median().round(2)}\t{df[3].mean().round(2)}\t{df[3].median().round(2)}')
    EOF

    python3 ./summary_stats.py \
      < ~{sample_id}.read_length_and_quality.tsv \
      > stats.txt

    cut -f1 stats.txt > num_reads.txt || echo "0" > num_reads.txt
    cut -f2 stats.txt > read_length_mean.txt || echo "0.00" > read_length_mean.txt
    cut -f3 stats.txt > read_length_median.txt || echo "0.00" > read_length_median.txt
    cut -f4 stats.txt > read_quality_mean.txt || echo "0.00" > read_quality_mean.txt
    cut -f5 stats.txt > read_quality_median.txt || echo "0.00" > read_quality_median.txt

    gzip ~{sample_id}.read_length_and_quality.tsv
  >>>

  output {
    File   read_length_and_quality  = "~{sample_id}.read_length_and_quality.tsv.gz"
    File   read_length_histogram    = "~{sample_id}.read_length_histogram.tsv"
    File   read_quality_histogram   = "~{sample_id}.read_quality_histogram.tsv"
    File   read_length_plot         = "~{sample_id}.read_length_histogram.png"
    File   read_quality_plot        = "~{sample_id}.read_quality_histogram.png"
    String stat_num_reads           = read_string("num_reads.txt")
    String stat_read_length_mean    = read_string("read_length_mean.txt")
    String stat_read_length_median  = read_string("read_length_median.txt")
    String stat_read_quality_mean   = read_string("read_quality_mean.txt")
    String stat_read_quality_median = read_string("read_quality_median.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}