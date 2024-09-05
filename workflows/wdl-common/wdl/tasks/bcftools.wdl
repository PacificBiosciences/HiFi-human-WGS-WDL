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

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcf, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools stats \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --apply-filters .,PASS \
      --samples ~{sample_id} \
      ~{"--fasta-ref " + ref_fasta} \
      ~{vcf} \
    > ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt

    # pull some top level stats
    grep -w '^SN' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | grep 'number of SNPs:' | cut -f4 > snv_count.txt
    grep -w '^SN' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | grep 'number of indels:' | cut -f4 > indel_count.txt
    grep -w '^TSTV' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | cut -f5 > tstv_ratio.txt
    nHets=$(grep -w '^PSC' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | cut -f6)
    nNonRefHom=$(grep -w '^PSC' ~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt | cut -f5)
    printf %.2f "$((10**2 * nHets / nNonRefHom))e-2" > hethom_ratio.txt  # hack for low precision float without bc

    bcftools roh \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --AF-dflt 0.4 \
      ~{vcf} \
      > ~{sample_id}.~{ref_name}.bcftools_roh.out

    # convert the roh output to a bed file with no filtering
    cat << EOF > roh_bed.py
    with open('~{sample_id}.~{ref_name}.bcftools_roh.out', 'r') as f:
      lines = f.readlines()
      print("#chr\tstart\tend\tqual")
      for line in lines:
        vals = line.strip().split('\t')
        if vals[0] == "RG":
          print('\t'.join([vals[2], vals[3], vals[4], vals[7]]))
    EOF

    python3 ./roh_bed.py > ~{sample_id}.~{ref_name}.roh.bed
  >>>

  output {
    File  stats              = "~{sample_id}.~{ref_name}.small_variants.vcf.stats.txt"
    File  roh_out            = "~{sample_id}.~{ref_name}.bcftools_roh.out"
    File  roh_bed            = "~{sample_id}.~{ref_name}.roh.bed"
    String stat_SNV_count    = read_string("snv_count.txt")
    String stat_INDEL_count  = read_string("indel_count.txt")
    String stat_TSTV_ratio   = read_string("tstv_ratio.txt")
    String stat_HETHOM_ratio = read_string("hethom_ratio.txt")
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
      ln -s "${input}" vcfs
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
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="DUP"' \
      "~{vcf}" | wc -l \
      > stat_DUP.txt || echo "0" > stat_DUP.txt
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="DEL"' \
      "~{vcf}" | wc -l \
      > stat_DEL.txt || echo "0" > stat_DEL.txt
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="INS"' \
      "~{vcf}" | wc -l \
      > stat_INS.txt || echo "0" > stat_INS.txt
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="INV"' \
      "~{vcf}" | wc -l \
      > stat_INV.txt || echo "0" > stat_INV.txt
    bcftools view -H -i 'FILTER="PASS" & SVTYPE="BND"' \
      "~{vcf}" | wc -l \
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