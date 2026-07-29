version 1.0

import "../structs.wdl"

task bcftools_stats_roh_small_variants {
  meta {
    description: "Run bcftools stats and bcftools roh on a small variant VCF"
    outputs: {
      stats: {
        description: "Small variant statistics"
      },
      roh_out: {
        description: "Regions of homozygosity"
      },
      roh_bed: {
        description: "Regions of homozygosity BED"
      },
      stat_SNV_count: {
        description: "Number of SNVs"
      },
      stat_INDEL_count: {
        description: "Number of INDELs"
      },
      stat_TSTV_ratio: {
        description: "Ts/Tv ratio"
      },
      stat_HETHOM_ratio: {
        description: "Het/Hom ratio for SNVs"
      },
      snv_distribution_plot: {
        description: "Distribution of SNVs by REF, ALT"
      },
      indel_distribution_plot: {
        description: "Distribution of indels by size"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    vcf: {
      description: "Small variant VCF"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_name: {
      description: "Reference name"
    }
    min_length: {
      description: "Minimum length for ROH"
    }
    min_qual: {
      description: "Minimum quality for ROH"
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
    File vcf
    File ref_fasta
    String ref_name
    Int min_length = 100000
    Int min_qual = 20
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    ln --symbolic --verbose "~{vcf}" "~{ref_fasta}" .

    bcftools norm \
      --fasta-ref "~{basename(ref_fasta)}" \
      --multiallelics - \
      "~{basename(vcf)}" 2>/dev/null \
    | bcftools view \
      --apply-filters .,PASS \
      --exclude 'GQ<20.0 || GT="ref" || GT="mis" || ALT="."' \
      --trim-alt-alleles \
      - \
    | bcftools stats \
      --samples "~{sample_id}" \
      --fasta-ref "~{basename(ref_fasta)}" \
      - \
    > "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt"

    # pull some top level stats
    grep -w '^SN' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" | grep 'number of SNPs:' | cut -f4 > snv_count.txt
    grep -w '^SN' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" | grep 'number of indels:' | cut -f4 > indel_count.txt
    grep -w '^TSTV' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" | cut -f5 > tstv_ratio.txt
    nHets=$(grep -w '^PSC' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" | cut -f6)
    nNonRefHom=$(grep -w '^PSC' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" | cut -f5)
    printf %.2f "$((10**2 * nHets / nNonRefHom))e-2" > hethom_ratio.txt  # hack for low precision float without bc

    # plot SNVs by REF and ALT, using the ST (substitution type) section bcftools stats already computed above
    cat << EOF > plot_snvs.py
    import sys, pandas as pd, seaborn as sns, matplotlib, numpy as np
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df = pd.read_csv(sys.stdin, sep='\t')
    df[['REF', 'ALT']] = df['type'].str.split('>', expand=True)
    df = pd.pivot(df, index='ALT', columns='REF', values='count')
    sns.set_style('dark')
    mask = np.identity(df.shape[0], dtype=bool)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(df, mask=mask, cmap='coolwarm', annot=True, fmt=',', annot_kws=dict(fontsize='large'), linewidth=.5, ax=ax)
    plt.xlabel('REF base', fontsize='large')
    plt.ylabel('ALT base', fontsize='large')
    plt.xticks(fontsize='large', rotation=0)
    plt.yticks(fontsize='large', rotation=0)
    plt.title('~{sample_id}.~{ref_name}\nDeepVariant SNV distribution', fontsize='large')
    fig.tight_layout()
    plt.savefig('~{sample_id}.~{ref_name}.small_variants.snv_distribution.png'); plt.close();
    EOF

    grep -w '^ST' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" \
    | cut -f3,4 \
    | awk -v OFS=$'\t' 'BEGIN {print "type", "count"} {print $1, $2}' \
    | python3 ./plot_snvs.py

    # plot indels by size, using the IDD (indel distribution) section bcftools stats already computed above
    cat << EOF > plot_indels.py
    import sys, pandas as pd, seaborn as sns, matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df = pd.read_csv(sys.stdin, sep='\t')
    def size_filter(df, col, min_len, max_len):
      return df[(df[col].abs() >= min_len) & (df[col].abs() < max_len)]
    def plot_hist(ax, df, min_len, max_len, logy=False, xlabel=True):
      g = sns.histplot(size_filter(df, 'length', min_len, max_len), x='length', weights='count', binwidth=1, binrange=(-max_len - 0.5, max_len + 0.5), ax=ax)
      g.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))
      g.set_xlim(-max_len, max_len)
      g.set_title(f'±[{min_len},{max_len}) bp')
      if not xlabel:
        g.set_xlabel('')
        g.tick_params(labelbottom=False)
      if logy:
        g.set_yscale('log')
        g.set_ylabel('log10(Count)')
    sns.set_style('darkgrid')
    fig, axs = plt.subplots(2, 1, figsize=(8,6))
    plot_hist(axs[0], df, 1, 50, xlabel=False)
    plot_hist(axs[1], df, 1, 50, logy=True)
    plt.suptitle('~{sample_id}.~{ref_name}\nDeepVariant indel distribution, ±[1,50) bp')
    fig.tight_layout()
    plt.savefig('~{sample_id}.~{ref_name}.small_variants.indel_distribution.png'); plt.close();
    EOF

    grep -w '^IDD' "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt" \
    | cut -f3,4 \
    | awk -v OFS=$'\t' 'BEGIN {print "length", "count"} {print $1, $2}' \
    | python3 ./plot_indels.py

    # find runs of homozygosity
    bcftools roh \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --AF-dflt 0.4 \
      "~{basename(vcf)}" \
    > "~{sample_id}.~{ref_name}.bcftools_roh.out"

    # convert the roh output to a bed file, filtering for length and quality
    cat << EOF > roh_bed.py
    with open('~{sample_id}.~{ref_name}.bcftools_roh.out', 'r') as f:
      print("#chr\tstart\tend\tqual")
      for line in f:
        if line.startswith("RG"):
          # RG [2]Sample [3]Chromosome [4]Start [5]End [6]Length (bp) [7]Number of markers [8]Quality (average fwd-bwd phred score)
          _, _, chrom, start, end, length, _, score = line.strip().split('\t')
          if int(length) >= ~{min_length} and float(score) >= ~{min_qual}:
            print('\t'.join([chrom, start, end, score]))
    EOF

    python3 ./roh_bed.py > "~{sample_id}.~{ref_name}.bcftools_roh.bed"

    gzip "~{sample_id}.~{ref_name}.bcftools_roh.out"
  >>>

  output {
    File stats = "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt"
    File roh_out = "~{sample_id}.~{ref_name}.bcftools_roh.out.gz"
    File roh_bed = "~{sample_id}.~{ref_name}.bcftools_roh.bed"
    String stat_SNV_count = read_string("snv_count.txt")
    String stat_INDEL_count = read_string("indel_count.txt")
    String stat_TSTV_ratio = read_string("tstv_ratio.txt")
    String stat_HETHOM_ratio = read_string("hethom_ratio.txt")
    File snv_distribution_plot = "~{sample_id}.~{ref_name}.small_variants.snv_distribution.png"
    File indel_distribution_plot = "~{sample_id}.~{ref_name}.small_variants.indel_distribution.png"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

task split_vcf_by_sample {
  meta {
    description: "Split a multi-sample VCF by sample"
    outputs: {
      split_vcfs: {
        description: "Split VCFs"
      },
      split_vcf_indices: {
        description: "Split VCF indices"
      }
    }
  }

  parameter_meta {
    sample_ids: {
      description: "Sample IDs"
    }
    vcf: {
      description: "VCF"
    }
    vcf_index: {
      description: "VCF index"
    }
    split_vcf_names: {
      description: "Split VCF names"
    }
    split_vcf_index_names: {
      description: "Split VCF index names"
    }
    exclude_uncalled: {
      description: "Exclude uncalled genotypes (default: true)"
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
    File vcf
    File vcf_index
    Array[String] split_vcf_names
    Array[String] split_vcf_index_names
    Boolean exclude_uncalled = true
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  String vcf_basename = basename(vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{vcf}" "~{vcf_index}" .

    bcftools --version

    for sample_id in ~{sep=" " sample_ids}; do
      # Extract sample, optionally excluding uncalled genotypes
        bcftools view \
          ~{if threads > 1
            then "--threads '" + (threads - 1) + "'"
            else ""
          } \
          --samples "${sample_id}" \
          ~{true="--exclude-uncalled" false="" exclude_uncalled} \
          --output-type z \
          --output "${sample_id}.~{vcf_basename}.vcf.gz" \
          --write-index=tbi \
          "~{basename(vcf)}"
    done
  >>>

  output {
    Array[File] split_vcfs = split_vcf_names
    Array[File] split_vcf_indices = split_vcf_index_names
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

task bcftools_merge {
  meta {
    description: "Merge multiple sample VCFs into a single joint VCF"
    outputs: {
      merged_vcf: {
        description: "Merged VCF"
      },
      merged_vcf_index: {
        description: "Merged VCF index"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "VCFs"
    }
    vcf_indices: {
      description: "VCF indices"
    }
    out_prefix: {
      description: "Output VCF name"
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
    Array[File] vcfs
    Array[File] vcf_indices
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    VCFS=()

    for i in ~{sep=" " vcfs}; do
      ln --symbolic --verbose "${i}" .
      # shellcheck disable=SC2086
      VCFS+=("$(basename ${i})")
    done
    for i in ~{sep=" " vcf_indices}; do
      ln --symbolic --verbose "${i}" .
    done

    bcftools --version

    # shellcheck disable=SC2068
    bcftools merge \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --output-type z \
      --output "~{out_prefix}.vcf.gz" \
      --write-index=tbi \
      ${VCFS[@]}
  >>>

  output {
    File merged_vcf = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

task sv_stats {
  meta {
    description: "Collect statistics on structural variants in a VCF"
    outputs: {
      stat_sv_DUP_count: {
        description: "Number of DUP structural variants"
      },
      stat_sv_DEL_count: {
        description: "Number of DEL structural variants"
      },
      stat_sv_INS_count: {
        description: "Number of INS structural variants"
      },
      stat_sv_INV_count: {
        description: "Number of INV structural variants"
      },
      stat_sv_BND_count: {
        description: "Number of BND structural variants"
      },
      stat_sv_SWAP_count: {
        description: "Number of structural variant sequence swap events"
      },
      sv_stats_plot: {
        description: "Structural variant size distribution plot"
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
    vcf: {
      description: "VCF"
    }
    min_length: {
      description: "Minimum length"
    }
    max_scar_length: {
      description: "Maximum scar length for INS/DEL to be considered simple"
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
    File vcf
    Int min_length = 50
    Int max_scar_length = 10
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    ln --symbolic --verbose "~{vcf}" .

    # Count the number of variants of each type. For DUP/DEL/INS/INV,
    # extract each variant's SVLEN via `bcftools query` instead of just
    # counting lines -- the same query then serves both the count stat
    # and the size-distribution plot below, with no extra VCF pass.
    bcftools query \
      -f '%INFO/SVLEN\n' \
      --include '(GT!="ref" & GT!="./." & GT!=".") & FILTER="PASS" & ABS(SVLEN)>=~{min_length} & SVTYPE="DUP"' \
      "~{basename(vcf)}" \
    > length_DUP.txt || : > length_DUP.txt
    wc --lines < length_DUP.txt > stat_DUP.txt
    bcftools query \
      -f '%INFO/SVLEN\n' \
      --include 'GT="alt" & FILTER="PASS" & ABS(SVLEN)>=~{min_length} & SVTYPE="DEL" & STRLEN(ALT[0])<=~{max_scar_length}' \
      "~{basename(vcf)}" \
    > length_DEL.txt || : > length_DEL.txt
    wc --lines < length_DEL.txt > stat_DEL.txt
    bcftools query \
      -f '%INFO/SVLEN\n' \
      --include 'GT="alt" & FILTER="PASS" & ABS(SVLEN)>=~{min_length} & SVTYPE="INS" & STRLEN(REF)<=~{max_scar_length}' \
      "~{basename(vcf)}" \
    > length_INS.txt || : > length_INS.txt
    wc --lines < length_INS.txt > stat_INS.txt
    bcftools view \
      --no-header \
      --include 'GT="alt" & FILTER="PASS" & ABS(SVLEN)>=~{min_length} & (SVTYPE="INS" | SVTYPE="DEL") & STRLEN(REF)>~{max_scar_length} & STRLEN(ALT[0])>~{max_scar_length}' \
      "~{basename(vcf)}" \
    | wc --lines \
    > stat_SWAP.txt || echo "0" > stat_SWAP.txt
    bcftools query \
      -f '%INFO/SVLEN\n' \
      --include 'GT="alt" & FILTER="PASS" & SVTYPE="INV"' \
      "~{basename(vcf)}" \
    > length_INV.txt || : > length_INV.txt
    wc --lines < length_INV.txt > stat_INV.txt
    bcftools view \
      --no-header \
      --include 'GT="alt" & FILTER="PASS" & SVTYPE="BND"' \
      "~{basename(vcf)}" \
    | wc --lines \
    > stat_BND.txt || echo "0" > stat_BND.txt

    # plot size distributions for INS/DEL (combined), DUP, and INV
    cat << EOF > plot_sv_stats.py
    import matplotlib, matplotlib.pyplot as plt
    matplotlib.use('Agg')
    import seaborn as sns
    from matplotlib.ticker import MaxNLocator, FuncFormatter
    sns.set_style('darkgrid')
    def fmt_bp(x, pos):
      sign = '-' if x < 0 else ''
      ax_val = abs(x)
      if ax_val >= 1000000:
        return f'{sign}{ax_val/1000000:g}M'
      if ax_val >= 1000:
        return f'{sign}{ax_val/1000:g}k'
      return f'{sign}{ax_val:g}'
    def read_lengths(path):
      with open(path) as f:
        values = []
        for line in f.read().split():
          for tok in line.split(','):
            try:
              values.append(int(tok))
            except ValueError:
              pass
        return values
    N_BINS_TOTAL = 50
    CAP = 2000000
    def tier_edges(lo, hi):
      # bin width derived as if this tier started at 0, then the bins that
      # would fall in [0, lo) are simply not displayed -- keeps bin width
      # (and bar count) consistent across tiers instead of shrinking bins
      # to fit each tier's own range
      w = hi // N_BINS_TOTAL
      return list(range(lo, hi, w)) + [hi]
    def plot_del_tier(ax, dele, lo, hi, ylabel=False, legend=False):
      # Plotted on absolute length (not signed), same tier order/direction
      # as the INS/DUP row directly below -- each figure column then shares
      # an identical x-axis across all three rows, rather than mirroring
      # DEL's sign against INS/DUP.
      cap = hi is None
      hi_eff = CAP if cap else hi
      edges = tier_edges(lo, hi_eff)
      vals = [min(-x, CAP) for x in dele if -x >= lo] if cap else [-x for x in dele if lo <= -x < hi]
      ax.hist(vals, bins=edges, color='C1', label='DEL')
      ax.yaxis.set_major_locator(MaxNLocator(integer=True))
      ax.set_xlabel('SV length (bp)')
      ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
      if ylabel:
        ax.set_ylabel('Count')
      if legend:
        ax.legend(fontsize='small')
      ax.set_title(f'>={lo:,} bp (capped {CAP:,})' if cap else f'[{lo:,}, {hi:,}) bp', fontsize='small')
    def plot_insdup_tier(ax, ins, dup, lo, hi, legend=False):
      # INS and DUP stacked together per tier (both are positive-length
      # events with no separate row of their own), ordered smallest-first
      # (leftmost, nearest zero) growing outward -- mirrors the DEL tiers
      cap = hi is None
      hi_eff = CAP if cap else hi
      edges = tier_edges(lo, hi_eff)
      if cap:
        ins_v = [min(x, CAP) for x in ins if x >= lo]
        dup_v = [min(x, CAP) for x in dup if x >= lo]
      else:
        ins_v = [x for x in ins if lo <= x < hi]
        dup_v = [x for x in dup if lo <= x < hi]
      ax.hist([ins_v, dup_v], bins=edges, stacked=True, color=['C0', 'C2'], label=['INS', 'DUP'])
      ax.yaxis.set_major_locator(MaxNLocator(integer=True))
      ax.set_xlabel('SV length (bp)')
      ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
      if legend:
        ax.legend(fontsize='small')
      ax.set_title(f'>={lo:,} bp (capped {CAP:,})' if cap else f'[{lo:,}, {hi:,}) bp', fontsize='small')
    def plot_inv_tier(ax, inv, lo, hi, legend=False):
      # single-series tiers, same layout/style as the INS/DUP row (smallest
      # nearest zero growing outward). INV isn't length-filtered elsewhere
      # in this task (unlike DUP/DEL/INS), so its own smallest tier starts
      # at 0, not min_length, to keep showing every INV regardless of size.
      cap = hi is None
      hi_eff = CAP if cap else hi
      edges = tier_edges(lo, hi_eff)
      vals = [min(x, CAP) for x in inv if x >= lo] if cap else [x for x in inv if lo <= x < hi]
      ax.hist(vals, bins=edges, color='C3', label='INV')
      ax.yaxis.set_major_locator(MaxNLocator(integer=True))
      ax.set_xlabel('SV length (bp)')
      ax.xaxis.set_major_formatter(FuncFormatter(fmt_bp))
      if legend:
        ax.legend(fontsize='small')
      ax.set_title(f'>={lo:,} bp (capped {CAP:,})' if cap else f'[{lo:,}, {hi:,}) bp', fontsize='small')
    ins = read_lengths('length_INS.txt')
    dele = read_lengths('length_DEL.txt')
    # SVLEN sign convention varies by caller/VCF version (e.g. sawfish 2.2.1 emits
    # positive SVLEN for DEL); force negative here so DEL always plots left of 0.
    dele = [-abs(x) for x in dele]
    dup = read_lengths('length_DUP.txt')
    inv = read_lengths('length_INV.txt')
    fig, axd = plt.subplot_mosaic(
      [['d1', 'd2', 'd3', 'd4'],
       ['i1', 'i2', 'i3', 'i4'],
       ['v1', 'v2', 'v3', 'v4']],
      figsize=(13, 12),
    )
    plot_del_tier(axd['d1'], dele, ~{min_length}, 1000, ylabel=True, legend=True)
    plot_del_tier(axd['d2'], dele, 1000, 10000)
    plot_del_tier(axd['d3'], dele, 10000, 100000)
    plot_del_tier(axd['d4'], dele, 100000, None)
    plot_insdup_tier(axd['i1'], ins, dup, ~{min_length}, 1000, legend=True)
    plot_insdup_tier(axd['i2'], ins, dup, 1000, 10000)
    plot_insdup_tier(axd['i3'], ins, dup, 10000, 100000)
    plot_insdup_tier(axd['i4'], ins, dup, 100000, None)
    plot_inv_tier(axd['v1'], inv, 0, 1000, legend=True)
    plot_inv_tier(axd['v2'], inv, 1000, 10000)
    plot_inv_tier(axd['v3'], inv, 10000, 100000)
    plot_inv_tier(axd['v4'], inv, 100000, None)
    # DEL and INS/DUP share a column's tier, so give each pair the same
    # y-limit -- bar heights become directly comparable, not just the
    # x-axis. INV keeps its own independent y-limits.
    for d_key, i_key in [('d1', 'i1'), ('d2', 'i2'), ('d3', 'i3'), ('d4', 'i4')]:
      shared_max = max(axd[d_key].get_ylim()[1], axd[i_key].get_ylim()[1])
      axd[d_key].set_ylim(0, shared_max)
      axd[i_key].set_ylim(0, shared_max)
    fig.suptitle('~{sample_id}.~{ref_name}\nStructural variant size distribution')
    fig.tight_layout()
    plt.savefig('~{sample_id}.~{ref_name}.sv_stats.png'); plt.close()
    EOF

    python3 ./plot_sv_stats.py
  >>>

  output {
    String stat_sv_DUP_count = read_string("stat_DUP.txt")
    String stat_sv_DEL_count = read_string("stat_DEL.txt")
    String stat_sv_INS_count = read_string("stat_INS.txt")
    String stat_sv_INV_count = read_string("stat_INV.txt")
    String stat_sv_BND_count = read_string("stat_BND.txt")
    String stat_sv_SWAP_count = read_string("stat_SWAP.txt")
    File sv_stats_plot = "~{sample_id}.~{ref_name}.sv_stats.png"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

task bcftools_rename {
  meta {
    description: "Renames sample IDs in a VCF using bcftools reheader"
    outputs: {
      renamed_vcf: {
        description: "Renamed VCF"
      },
      renamed_vcf_index: {
        description: "Renamed VCF index"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "Input VCF"
    }
    sample_id_lookup: {
      description: "Sample ID lookup table, as an array of [old_id, new_id] arrays"
    }
    out_prefix: {
      description: "Output VCF prefix"
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
    File vcf
    Array[Array[String]] sample_id_lookup
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    ln --symbolic --verbose "~{vcf}" .

    bcftools reheader \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --samples "~{write_tsv(sample_id_lookup)}" \
      --output "~{out_prefix}.vcf.gz" \
      "~{basename(vcf)}"
    bcftools index --tbi \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
    "~{out_prefix}.vcf.gz"
  >>>

  output {
    File renamed_vcf = "~{out_prefix}.vcf.gz"
    File renamed_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # pb_wdl_base:build4
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

