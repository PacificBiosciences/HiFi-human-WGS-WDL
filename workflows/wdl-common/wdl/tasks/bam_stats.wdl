version 1.0

import "../structs.wdl"

task bam_stats {
  meta {
    description: "Get read, mapping, and alignment statistics from a BAM file."
    outputs: {
      bam_statistics: {
        description: "BAM statistics"
      },
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
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String ref_name
    File bam
    File bam_index
    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(bam, "GB") + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{bam}" .
    ln --symbolic --verbose "~{bam_index}" .

    # pull the count and percent of mapped reads
    samtools flagstat \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --output-fmt json \
      "~{basename(bam)}" \
    | jq '.["QC-passed reads"]["primary mapped %", "primary mapped"]' \
    > mapping_stats.txt
    sed -n '1p' mapping_stats.txt > mapped_percent.txt
    sed -n '2p' mapping_stats.txt > mapped_read_count.txt

    cat << EOF > extract_stats.py
    import math, pysam
    MAX_QV = 60
    save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
    with pysam.AlignmentFile('~{basename(bam)}', check_sq=False) as bamin:
      pysam.set_verbosity(save)  # restore warnings
      for r in bamin:
        if r.is_unmapped:
          status = 'unmapped'
        elif r.is_supplementary:
          status = 'supp'
        else:
          status = 'prim'
        if math.isnan(rq := r.get_tag('rq') if r.has_tag('rq') else math.nan):
          readqv = math.nan
        else:
          readqv = MAX_QV if (errorrate := 1.0 - rq) == 0 else math.floor(-10 * math.log10(errorrate))
        print("\t".join([
          r.query_name.split('/')[0],
          r.query_name,
          str(len(r.query_sequence)),
          str(readqv),
          status,
          str(math.nan if r.is_unmapped else r.mapping_quality),
          f"{math.nan if r.is_unmapped else r.get_tag('mg'):.4f}"
        ]))
    EOF

    cat << EOF > plot_stats.py
    import pandas as pd, seaborn as sns, matplotlib, matplotlib.pyplot as plt
    matplotlib.use('Agg')
    df = pd.read_csv(
      "~{sample_id}.~{ref_name}.bam_statistics.tsv.gz", sep='\t', header=None, compression='gzip',
      names=['movie', 'read_name', 'read_length', 'read_quality', 'alignment_type', 'mapq', 'mg'],
      usecols=['movie', 'read_length', 'read_quality', 'alignment_type', 'mapq', 'mg'],
      dtype={
        'movie': 'category', 'read_length': pd.Int32Dtype(), 'read_quality': pd.Int16Dtype(),
        'alignment_type': 'category', 'mapq': pd.Int8Dtype(), 'mg': pd.Float32Dtype(),
      },
    )
    prim = df['alignment_type'] == 'prim'      # primary alignments
    aligned = df['alignment_type'].isin(['prim', 'supp'])  # all aligned reads
    non_supp = df['alignment_type'] != 'supp'  # unique reads (primary and unmapped)
    has_quality = df['read_quality'].notna()   # reads with quality scores
    # calculate length n50
    read_lengths = df[non_supp].read_length.sort_values(ascending=False).reset_index(drop=True)
    cumsum = read_lengths.cumsum()
    read_length_n50 = read_lengths[cumsum >= cumsum.iloc[-1] / 2].iloc[0]
    # print summary stats to stdout
    print('\t'.join([
      str(len(df[non_supp])),
      str(len(df[prim])),
      str(round(len(df[prim]) / len(df[non_supp]) * 100, 2)),
      str(df[non_supp].read_length.mean().round(2)),
      str(df[non_supp].read_length.median().round(2)),
      str(read_length_n50),
      str(df[non_supp].read_quality.mean().round(2) if any(has_quality) else pd.NA),
      str(df[non_supp].read_quality.median().round(2) if any(has_quality) else pd.NA),
      str(df[aligned].mg.mean().round(2)),
      str(df[aligned].mg.median().round(2))
    ]))
    sns.set_theme(style='darkgrid')
    count_formatter = plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}')
    # plot read length histogram
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(data=df[non_supp], x='read_length', hue='movie', bins=range(0, 40000, 1000), multiple='stack', ax=ax)
    ax.set_xlim(0, 40000); ax.set_xlabel('read length (bp)'); ax.set_ylabel('read count');
    ax.yaxis.set_major_formatter(count_formatter);
    ax.set_title('~{sample_id}\nRead length histogram'); fig.tight_layout();
    plt.savefig('~{sample_id}.read_length_histogram.png'); plt.close();
    # plot read quality histogram
    if len(df[has_quality & non_supp]) > 0:
      fig, ax = plt.subplots(figsize=(8, 6))
      sns.histplot(data=df[has_quality & non_supp], x='read_quality', hue='movie', bins=range(0, 61), multiple='stack', ax=ax)
      ax.set_xlabel('Phred-scaled read quality'); ax.set_xlim(0, 60); ax.set_ylabel('read count');
      ax.yaxis.set_major_formatter(count_formatter);
      ax.set_title('~{sample_id}\nRead quality histogram'); fig.tight_layout();
      plt.savefig('~{sample_id}.read_quality_histogram.png'); plt.close();
    # plot mapq histograms
    fig, axs = plt.subplots(2, 1, figsize=(8, 6))
    sns.histplot(df, x='mapq', hue='alignment_type', element='step', kde=False, bins=60, ax=axs[0])
    sns.histplot(df, x='mapq', hue='alignment_type', element='step', kde=False, bins=60, ax=axs[1], legend=False)
    sns.move_legend(axs[0], 'upper left');
    axs[0].axes.set_xlim(0, 60); axs[0].axes.set_xlabel(''); axs[0].axes.set_xticklabels([]);
    axs[1].axes.set_xlim(0, 60); axs[1].axes.set_yscale('log'); axs[1].axes.set_ylabel('log10(Count)');
    fig.suptitle('~{sample_id}.~{ref_name}\nMAPQ distribution'); plt.tight_layout();
    plt.savefig('~{sample_id}.~{ref_name}.mapq_distribution.png'); plt.close();
    # plot mg histograms
    fig, axs = plt.subplots(2, 1, figsize=(8, 6))
    sns.histplot(df, x='mg', hue='alignment_type', element='step', kde=False, bins=60, ax=axs[0])
    sns.histplot(df, x='mg', hue='alignment_type', element='step', kde=False, bins=60, ax=axs[1], legend=False)
    sns.move_legend(axs[0], 'upper left');
    axs[0].set_xlim(70, 100); axs[0].axes.set_xlabel(''); axs[0].axes.set_xticklabels([]);
    axs[1].set_xlim(70, 100); axs[1].axes.set_yscale('log'); axs[1].axes.set_ylabel('log10(Count)');
    plt.suptitle('~{sample_id}.~{ref_name}\nGap-compressed sequence identity distribution'); fig.tight_layout();
    plt.savefig('~{sample_id}.~{ref_name}.mg_distribution.png'); plt.close();
    EOF

    python3 extract_stats.py \
    | gzip --stdout > "~{sample_id}.~{ref_name}.bam_statistics.tsv.gz"

    python3 ./plot_stats.py > stats.txt

    cut -f1 stats.txt > read_count.txt
    cut -f2 stats.txt > mapped_read_count.txt
    cut -f3 stats.txt > mapped_read_percent.txt
    cut -f4 stats.txt > read_length_mean.txt
    cut -f5 stats.txt > read_length_median.txt
    cut -f6 stats.txt > read_length_n50.txt
    cut -f7 stats.txt > read_quality_mean.txt
    cut -f8 stats.txt > read_quality_median.txt
    cut -f9 stats.txt > mg_mean.txt
    cut -f10 stats.txt > mg_median.txt
  >>>

  output {
    File bam_statistics = "~{sample_id}.~{ref_name}.bam_statistics.tsv.gz"
    File read_length_plot = "~{sample_id}.read_length_histogram.png"
    File? read_quality_plot = "~{sample_id}.read_quality_histogram.png"
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
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
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
