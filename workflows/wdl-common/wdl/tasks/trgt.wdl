version 1.0

import "../structs.wdl"

task trgt {
  meta {
    description: "Genotype tandem repeats from aligned reads using TRGT"
    outputs: {
      bam: {
        description: "Aligned TRGT spanning reads"
      },
      bam_index: {
        description: "Index for aligned TRGT spanning reads"
      },
      vcf: {
        description: "Phased TRGT VCF"
      },
      vcf_index: {
        description: "Index for phased TRGT VCF"
      },
      dropouts: {
        description: "TRGT regions with coverage dropouts"
      },
      stat_genotyped_count: {
        description: "Number of sites genotyped by TRGT"
      },
      stat_uncalled_count: {
        description: "Number of sites ungenotyped by TRGT"
      },
      msg: {
        description: "Array of messages"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    sex: {
      description: "Sample sex",
      choices: [
        "MALE",
        "FEMALE"
      ]
    }
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Aligned BAM index"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    trgt_bed: {
      description: "TRGT tandem repeat catalog BED"
    }
    expected_male_bed: {
      description: "Expected ploidy BED for sample with XY karyotype"
    }
    expected_female_bed: {
      description: "Expected ploidy BED for sample with XX karyotype"
    }
    out_prefix: {
      description: "Output prefix"
    }
    max_depth: {
      description: "Maximum locus depth"
    }
    min_read_quality: {
      description: "Minimum read quality"
    }
    haplotype_coverage_threshold: {
      description: "Haplotype coverage threshold for dropout analysis"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String? sex
    File aligned_bam
    File aligned_bam_index
    File ref_fasta
    File ref_index
    File trgt_bed
    File expected_male_bed
    File expected_female_bed
    String out_prefix
    Int max_depth = 50
    Float min_read_quality = 0.98
    Int haplotype_coverage_threshold = 2
    RuntimeAttributes runtime_attributes
  }

  Int threads = 32
  Int mem_gb = 64
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  Int samtools_sort_threads = 8

  String karyotype = if select_first([
    sex,
    "FEMALE"
  ]) == "MALE"
    then "XY"
    else "XX"
  File expected_bed = if select_first([
    sex,
    "FEMALE"
  ]) == "MALE"
    then expected_male_bed
    else expected_female_bed

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{aligned_bam}" .
    ln --symbolic --verbose "~{aligned_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .
    ln --symbolic --verbose "~{trgt_bed}" .
    ln --symbolic --verbose "~{expected_bed}" .

    if [ "~{defined(sex)}" != "true" ]; then
      echo "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for TRGT." >> messages.txt
    fi

    trgt --version

    trgt genotype \
      --threads ~{threads} \
      --karyotype "~{karyotype}" \
      --genome "~{basename(ref_fasta)}" \
      --repeats "~{basename(trgt_bed)}" \
      --reads "~{basename(aligned_bam)}" \
      --max-depth ~{max_depth} \
      --min-read-quality=~{min_read_quality} \
      --output-prefix "~{out_prefix}.trgt"

    bcftools --version

    bcftools sort \
      --output-type z \
      --output "~{out_prefix}.trgt.sorted.vcf.gz" \
      "~{out_prefix}.trgt.vcf.gz"

    bcftools index \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --tbi \
      "~{out_prefix}.trgt.sorted.vcf.gz"

    samtools --version

    samtools sort \
      --threads ~{samtools_sort_threads} \
      -m 800M \
      -o "~{out_prefix}.trgt.spanning.sorted.bam" \
      "~{out_prefix}.trgt.spanning.bam"
    samtools index \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      "~{out_prefix}.trgt.spanning.sorted.bam"

    find_trgt_dropouts.py \
      --ploidybed "~{basename(expected_bed)}" \
      --coverage ~{haplotype_coverage_threshold} \
      "~{basename(trgt_bed)}" \
      "~{out_prefix}.trgt.spanning.sorted.bam" \
      > "~{out_prefix}.trgt.dropouts.txt"

    bcftools view --no-header --exclude-uncalled \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      "~{out_prefix}.trgt.sorted.vcf.gz" \
    | wc --lines > genotyped_count.txt || echo "0" > genotyped_count.txt

    bcftools view --no-header --uncalled \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      "~{out_prefix}.trgt.sorted.vcf.gz"   \
    | wc --lines > uncalled_count.txt || echo "0" > uncalled_count.txt
  >>>

  output {
    File bam = "~{out_prefix}.trgt.spanning.sorted.bam"
    File bam_index = "~{out_prefix}.trgt.spanning.sorted.bam.bai"
    File vcf = "~{out_prefix}.trgt.sorted.vcf.gz"
    File vcf_index = "~{out_prefix}.trgt.sorted.vcf.gz.tbi"
    File dropouts = "~{out_prefix}.trgt.dropouts.txt"
    String stat_genotyped_count = read_string("genotyped_count.txt")
    String stat_uncalled_count = read_string("uncalled_count.txt")
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:be0ed7c173d221bd84e360b2b056e2abbecadd07ed86ffd4883a5cecca7a1e57"  # 5.0.0_build3
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

task trgt_merge {
  meta {
    description: "Merge results multiple samples analyzed with TRGT"
    outputs: {
      merged_vcf: {
        description: "Merged TRGT VCF"
      },
      merged_vcf_index: {
        description: "Index for merged TRGT VCF"
      }
    }
  }

  parameter_meta {
    vcfs: {
      description: "TRGT VCFs"
    }
    vcf_indices: {
      description: "TRGT VCF indices"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    out_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices
    File ref_fasta
    File ref_index
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb = 8
  Int disk_size = ceil((size(vcfs, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    trgt --version

    # shellcheck disable=SC2086
    trgt merge \
      --vcf ~{sep=" " vcfs} \
      --genome "~{basename(ref_fasta)}" \
      --output "~{out_prefix}.vcf.gz" \
      --write-index \
      --no-index
  >>>

  output {
    File merged_vcf = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:be0ed7c173d221bd84e360b2b056e2abbecadd07ed86ffd4883a5cecca7a1e57"  # 5.0.0_build3
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
