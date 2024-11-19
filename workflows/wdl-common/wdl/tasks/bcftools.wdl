version 1.0

import "../structs.wdl"

task bcftools_stats_roh_small_variants {
  meta {
    description: "Run bcftools stats and bcftools roh on a small variant VCF."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    vcf: {
      name: "Small variant VCF"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_name: {
      name: "Reference name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    stats: {
      name: "Small variant VCF stats"
    }
    roh_out: {
      name: "Runs of homozygosity output"
    }
    roh_bed: {
      name: "Runs of homozygosity BED"
    }
    stat_SNV_count: {
      name: "SNV count"
    }
    stat_INDEL_count: {
      name: "INDEL count"
    }
    stat_TSTV_ratio: {
      name: "Ts/Tv ratio"
    }
    stat_HETHOM_ratio: {
      name: "SNV Het/Hom ratio"
    }
  }

  input {
    String sample_id

    File vcf

    File ref_fasta
    String ref_name

    Int min_length = 100000
    Int min_qual = 20

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcf, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools norm \
      --fasta-ref ~{ref_fasta} \
      --multiallelics - \
      ~{vcf} 2>/dev/null \
    | bcftools view \
      --apply-filters .,PASS \
      --exclude 'GQ<20.0 || GT="ref" || GT="mis" || ALT="."' \
      --trim-alt-alleles \
      - \
    | bcftools stats \
      --samples ~{sample_id} \
      ~{"--fasta-ref " + ref_fasta} \
      - \
    > ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt

    # pull some top level stats
    grep -w '^SN' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | grep 'number of SNPs:' | cut -f4 > snv_count.txt
    grep -w '^SN' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | grep 'number of indels:' | cut -f4 > indel_count.txt
    grep -w '^TSTV' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | cut -f5 > tstv_ratio.txt
    nHets=$(grep -w '^PSC' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | cut -f6)
    nNonRefHom=$(grep -w '^PSC' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | cut -f5)
    printf %.2f "$((10**2 * nHets / nNonRefHom))e-2" > hethom_ratio.txt  # hack for low precision float without bc

    # plot SNVs by REF and ALT
    cat << EOF > plot_snvs.py
    import sys, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, numpy as np
    BASES = ['A', 'C', 'G', 'T']
    df = pd.concat(
      [
        pd.read_csv(sys.stdin, sep='\t'),
        pd.DataFrame.from_dict({'REF': BASES, 'ALT': BASES, 'count': [0] * len(BASES)})
      ],
      ignore_index=True
    )
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
    plt.savefig('~{sample_id}.~{ref_name}.small_variants.snv_distribution.png')
    EOF

    # normalize VCF, filtering for PASS SNVs >= GQ20, group by REF and ALT, and plot
    bcftools norm \
      --fasta-ref ~{ref_fasta} \
      --multiallelics - \
      ~{vcf} 2>/dev/null \
    | bcftools view \
      --apply-filters .,PASS \
      --include 'TYPE="snp" && GQ>=20.0 && GT="alt"' \
      --trim-alt-alleles \
      - \
    | bcftools query \
      --format '%REF\t%ALT\n' \
      - \
    | sort | uniq -c | sed 's/^\s*//;s/\s/\t/' \
    | awk -v OFS=$'\t' 'BEGIN {print "REF", "ALT", "count"} {print $2, $3, $1}' \
    | python3 ./plot_snvs.py -

    # plot indels by size
    cat << EOF > plot_indels.py
    import sys, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
    from numpy import abs
    df = pd.read_csv(sys.stdin, sep='\t')
    def size_filter(df, col, min, max):
      return df[(abs(df[col]) >= min) & (abs(df[col]) < max)]
    def plot_hist(ax, df, min, max, bins=100, logy=False, xlabel=True):
      g = sns.histplot(size_filter(df, 'length', min, max), x='length', weights='count', bins=bins, ax=ax)
      g.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))
      g.set_xlim(-max, max)
      g.set_title(f'±[{min},{max}) bp')
      if not xlabel:
        g.set_xlabel('')
        g.set_xticklabels([])
      if logy:
        g.set_yscale('log')
        g.set_ylabel('log10(Count)')
    sns.set_style('darkgrid')
    fig, axs = plt.subplots(2, 1, figsize=(8,6))
    plot_hist(axs[0], df, 1, 50, xlabel=False)
    plot_hist(axs[1], df, 1, 50, logy=True)
    plt.suptitle('~{sample_id}.~{ref_name}\nDeepVariant indel distribution, ±[1,50) bp')
    fig.tight_layout()
    plt.savefig('~{sample_id}.~{ref_name}.small_variants.indel_distribution.png')
    EOF

    # normalize VCF, filter for PASS indels >= GQ20, group by length, and plot
    bcftools norm \
      --fasta-ref ~{ref_fasta} \
      --multiallelics - \
      ~{vcf} 2>/dev/null \
    | bcftools view \
      --apply-filters .,PASS \
      --include 'TYPE="indel" && GQ>=20.0' \
      --trim-alt-alleles \
      - \
    | bcftools query \
      --format '%REF\t%ALT\n' \
      - \
    | awk -v OFS=$'\t' '{print length($2) - length($1)}' \
    | sort -n | uniq -c | sed 's/^\s*//' \
    | awk -v OFS=$'\t' 'BEGIN {print "length", "count"} {print $2, $1}' \
    | python3 ./plot_indels.py -

    # find runs of homozygosity
    bcftools roh \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --AF-dflt 0.4 \
      ~{vcf} \
    > ~{sample_id}.~{ref_name}.bcftools_roh.out

    # convert the roh output to a bed file, filtering for length and quality
    cat << EOF > roh_bed.py
    with open('~{sample_id}.~{ref_name}.bcftools_roh.out', 'r') as f:
      lines = f.readlines()
      print("#chr\tstart\tend\tqual")
      for line in lines:
        if line.startswith("RG"):
          # RG [2]Sample [3]Chromosome [4]Start [5]End [6]Length (bp) [7]Number of markers [8]Quality (average fwd-bwd phred score)
          _, _, chr, start, end, length, _, score = line.strip().split('\t')
          if int(length) >= ~{min_length} and float(score) >= ~{min_qual}:
            print('\t'.join([chr, start, end, score]))
    EOF

    python3 ./roh_bed.py > ~{sample_id}.~{ref_name}.bcftools_roh.bed

    gzip "~{sample_id}.~{ref_name}.bcftools_roh.out"
  >>>

  output {
    File   stats                   = "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt"
    File   roh_out                 = "~{sample_id}.~{ref_name}.bcftools_roh.out.gz"
    File   roh_bed                 = "~{sample_id}.~{ref_name}.bcftools_roh.bed"
    String stat_SNV_count          = read_string("snv_count.txt")
    String stat_INDEL_count        = read_string("indel_count.txt")
    String stat_TSTV_ratio         = read_string("tstv_ratio.txt")
    String stat_HETHOM_ratio       = read_string("hethom_ratio.txt")
    File   snv_distribution_plot   = "~{sample_id}.~{ref_name}.small_variants.snv_distribution.png"
    File   indel_distribution_plot = "~{sample_id}.~{ref_name}.small_variants.indel_distribution.png"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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

task concat_pbsv_vcf {
  meta {
    description: "Concatenate multiple PBSV VCFs into a single VCF."
  }

  parameter_meta {
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    out_prefix: {
      name: "Output VCF prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    concatenated_vcf: {
      name: "Concatenated VCF"
    }
    concatenated_vcf_index: {
      name: "Concatenated VCF index"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    mkdir vcfs
    while read -r input || [[ -n "${input}" ]]; do
      ln --verbose --symbolic "${input}" vcfs
    done < ~{write_lines(flatten([vcfs,vcf_indices]))}

    find vcfs -name "*.vcf.gz" > vcf.list

    bcftools --version

    bcftools concat \
      --allow-overlaps \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz \
      --file-list vcf.list
    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File concatenated_vcf       = "~{out_prefix}.vcf.gz"
    File concatenated_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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

task split_vcf_by_sample {
  meta {
    description: "Split a multi-sample VCF by sample."
  }

  parameter_meta {
    sample_ids: {
      name: "Sample IDs"
    }
    vcf: {
      name: "VCF"
    }
    vcf_index: {
      name: "VCF index"
    }
    split_vcf_names: {
      name: "Split VCF names"
    }
    split_vcf_index_names: {
      name: "Split VCF index names"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    split_vcfs: {
      name: "Split VCFs"
    }
    split_vcf_indices: {
      name: "Split VCF indices"
    }
  }

  input {
    Array[String] sample_ids
    File vcf
    File vcf_index

    Array[String] split_vcf_names
    Array[String] split_vcf_index_names

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  String vcf_basename = basename(vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    bcftools --version

    for sample_id in ~{sep=" " sample_ids}; do
      # Extract sample, keeping only passing variants and excluding uncalled genotypes
      bcftools view \
        ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
        --samples ${sample_id} \
        --exclude-uncalled \
        --output-type z \
        --output ${sample_id}.~{vcf_basename}.vcf.gz \
        ~{vcf}
      bcftools index --tbi \
        ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
        ${sample_id}.~{vcf_basename}.vcf.gz
      echo ${sample_id}.~{vcf_basename}.vcf.gz >> vcf.list
      echo ${sample_id}.~{vcf_basename}.vcf.gz.tbi >> index.list
    done
  >>>

  output {
    Array[File] split_vcfs        = split_vcf_names
    Array[File] split_vcf_indices = split_vcf_index_names
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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

task bcftools_merge {
  meta {
    description: "Merge multiple sample VCFs into a single joint VCF."
  }

  parameter_meta {
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    out_prefix: {
      name: "Output VCF name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_vcf: {
      name: "Merged VCF"
    }
    merged_vcf_index: {
      name: "Merged VCF index"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools merge \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz \
      ~{sep=" " vcfs}
    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File merged_vcf       = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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

task sv_stats {
  meta {
    description: "Collect statistics on structural variants in a VCF."
  }

  parameter_meta {
    vcf: {
      name: "VCF"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    stat_sv_DUP_count: {
      name: "Number of passing duplications greater than 49 bp"
    }
    stat_sv_DEL_count: {
      name: "Number of passing deletions greater than 49 bp"
    }
    stat_sv_INS_count: {
      name: "Number of passing insertions greater than 49 bp"
    }
    stat_sv_INV_count: {
      name: "Number of passing inversions greater than 49 bp"
    }
    stat_sv_BND_count: {
      name: "Number of breakends"
    }
  }

  input {
    File vcf

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcf, "GB") + 20)

  command <<<
    # Count the number of variants of each type
    bcftools view \
      --no-header \
      --include 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="DUP"' \
      "~{vcf}" \
    | wc --lines \
    > stat_DUP.txt || echo "0" > stat_DUP.txt
    bcftools view \
      --no-header \
      --include 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="DEL"' \
      "~{vcf}" \
    | wc --lines \
    > stat_DEL.txt || echo "0" > stat_DEL.txt
    bcftools view \
      --no-header \
      --include 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="INS"' \
      "~{vcf}" \
    | wc --lines \
    > stat_INS.txt || echo "0" > stat_INS.txt
    bcftools view \
      --no-header \
      --include 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="INV"' \
      "~{vcf}" \
    | wc --lines \
    > stat_INV.txt || echo "0" > stat_INV.txt
    bcftools view \
      --no-header \
      --include 'FILTER="PASS" & SVTYPE="BND"' \
      "~{vcf}" \
    | wc --lines \
    > stat_BND.txt || echo "0" > stat_BND.txt
  >>>

  output {
    String stat_sv_DUP_count = read_string("stat_DUP.txt")
    String stat_sv_DEL_count = read_string("stat_DEL.txt")
    String stat_sv_INS_count = read_string("stat_INS.txt")
    String stat_sv_INV_count = read_string("stat_INV.txt")
    String stat_sv_BND_count = read_string("stat_BND.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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