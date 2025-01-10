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

    cat << EOF > plot_length_and_quality.py
    import sys, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
    from numpy import nan
    sns.set_theme(style='darkgrid')
    df = pd.read_csv(sys.stdin, sep='\t', header=None,
      names=['movie', 'read_name', 'read_length', 'read_quality'],
      dtype={'movie': str, 'read_name': str, 'read_length': pd.Int32Dtype(), 'read_quality': pd.Int16Dtype()})
    len_mean = df["read_length"].mean().round(2)
    len_median = df["read_length"].median().round(2)
    qual_mean = df["read_quality"].mean().round(2) if not all(df["read_quality"].isna()) else nan
    qual_median = df["read_quality"].median().round(2) if not all(df["read_quality"].isna()) else nan
    print(f'{len(df)}\t{len_mean}\t{len_median}\t{qual_mean}\t{qual_median}')
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(data=df, x='read_length', hue='movie', bins=range(0, 40000, 1000), multiple='stack', ax=ax)
    ax.set_xlim(0,40000); ax.set_xlabel('read length (bp)');
    ax.set_ylabel('read count'); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'));
    ax.set_title('~{sample_id}\nRead length histogram'); fig.tight_layout();
    plt.savefig('~{sample_id}.read_length_histogram.png')
    df = df[df['read_quality'].notna()]
    if len(df) > 0:
      fig, ax = plt.subplots(figsize=(8, 6))
      sns.histplot(data=df, x='read_quality', hue='movie', bins=range(0, 61), multiple='stack', ax=ax)
      ax.set_xlabel('Phred-scaled read quality'); ax.set_xlim(0,60);
      ax.set_ylabel('read count'); ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'));
      ax.set_title('~{sample_id}\nRead quality histogram'); fig.tight_layout();
      plt.savefig('~{sample_id}.read_quality_histogram.png')
    EOF

    python3 ./plot_length_and_quality.py < ~{sample_id}.read_length_and_quality.tsv > stats.txt

    cut -f1 stats.txt > num_reads.txt
    cut -f2 stats.txt > read_length_mean.txt
    cut -f3 stats.txt > read_length_median.txt
    cut -f4 stats.txt > read_quality_mean.txt
    cut -f5 stats.txt > read_quality_median.txt

    gzip ~{sample_id}.read_length_and_quality.tsv
  >>>

  output {
    File   read_length_and_quality  = "~{sample_id}.read_length_and_quality.tsv.gz"
    File   read_length_plot         = "~{sample_id}.read_length_histogram.png"
    File?  read_quality_plot        = "~{sample_id}.read_quality_histogram.png"
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