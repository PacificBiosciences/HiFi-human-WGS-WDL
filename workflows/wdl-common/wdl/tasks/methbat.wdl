version 1.0

import "../structs.wdl"

task methbat_pileup {
  meta {
    description: "Generate methylation pileups from aligned reads"
    outputs: {
      cpg_pileup_bed: {
        description: "5mCpG pileup BED"
      },
      cpg_pileup_bed_index: {
        description: "Index for 5mCpG pileup BED"
      },
      hmcpg_pileup_bed: {
        description: "5hmC pileup BED"
      },
      hmcpg_pileup_bed_index: {
        description: "Index for 5hmC pileup BED"
      },
      ma_pileup_bed: {
        description: "6mA pileup BED"
      },
      ma_pileup_bed_index: {
        description: "Index for 6mA pileup BED"
      },
      summary: {
        description: "Summary JSON with counts of profiled sites and other metrics"
      },
      stat_hap1_cpg_count: {
        description: "Number of scored 5mCpGs in haplotype 1"
      },
      stat_hap2_cpg_count: {
        description: "Number of scored 5mCpGs in haplotype 2"
      },
      stat_combined_cpg_count: {
        description: "Number of scored 5mCpGs combined"
      },
      msg: {
        description: "Messages from MethBat execution, including errors if the command failed"
      }
    }
  }

  parameter_meta {
    haplotagged_bam: {
      description: "Aligned BAM"
    }
    haplotagged_bam_index: {
      description: "Aligned BAM index"
    }
    min_mapq: {
      description: "Minimum mapping quality to consider a read"
    }
    min_coverage: {
      description: "Minimum number of reads with a modification tag required to report a modification site"
    }
    edge_trimming_size: {
      description: "Number of bases to trim from the edges of the reads, for both coverage and methylation values"
    }
    phase_set_min_fraction: {
      description: "Minimum fraction of haplotagged reads a phase set must represent to be treated as dominant"
    }
    include_empty_sites: {
      description: "Enables the reporting of sites with no methylated reads for 6mA"
    }
    skip_5mC: {
      description: "Skip 5mC pileup generation"
    }
    skip_5hmC: {
      description: "Skip 5hmC pileup generation"
    }
    skip_6mA: {
      description: "Skip 6mA pileup generation"
    }
    out_prefix: {
      description: "Output prefix"
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
    File haplotagged_bam
    File haplotagged_bam_index
    Int min_mapq = 1
    Int min_coverage = 4
    Int edge_trimming_size = 20
    Float phase_set_min_fraction = 0.7
    Boolean include_empty_sites = false
    Boolean skip_5mC = false
    Boolean skip_5hmC = false
    Boolean skip_6mA = true
    String out_prefix
    Int threads = 8
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(haplotagged_bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    touch messages.txt

    ln --symbolic --verbose "~{haplotagged_bam}" "~{haplotagged_bam_index}" .

    methbat --version

    methbat pileup \
      --threads ~{threads} \
      --input-bam "~{basename(haplotagged_bam)}" \
      --min-mapq ~{min_mapq} \
      --min-coverage ~{min_coverage} \
      --edge-trimming-size ~{edge_trimming_size} \
      --phase-set-min-fraction ~{phase_set_min_fraction} \
      ~{if (include_empty_sites)
        then "--include-empty-sites"
        else ""
      } \
      ~{if (skip_5mC)
        then "--skip-5mC"
        else ""
      } \
      ~{if (skip_5hmC)
        then "--skip-5hmC"
        else ""
      } \
      ~{if (skip_6mA)
        then "--skip-6mA"
        else ""
      } \
      --output-prefix "~{out_prefix}" \
    || echo "MethBat pileup failed" >> messages.txt

    echo "0" > "~{out_prefix}.combined.bed.count"
    echo "0" > "~{out_prefix}.hap1.bed.count"
    echo "0" > "~{out_prefix}.hap2.bed.count"

    if [ -f "~{out_prefix}.5mC.bed.gz" ]; then
      zgrep -v '^#' "~{out_prefix}.5mC.bed.gz" | grep -c 'Total' > "~{out_prefix}.combined.bed.count" || true
      zgrep -v '^#' "~{out_prefix}.5mC.bed.gz" | grep -c 'hap1' > "~{out_prefix}.hap1.bed.count" || true
      zgrep -v '^#' "~{out_prefix}.5mC.bed.gz" | grep -c 'hap2' > "~{out_prefix}.hap2.bed.count" || true
    fi
  >>>

  output {
    File? cpg_pileup_bed = "~{out_prefix}.5mC.bed.gz"
    File? cpg_pileup_bed_index = "~{out_prefix}.5mC.bed.gz.tbi"
    File? hmcpg_pileup_bed = "~{out_prefix}.5hmC.bed.gz"
    File? hmcpg_pileup_bed_index = "~{out_prefix}.5hmC.bed.gz.tbi"
    File? ma_pileup_bed = "~{out_prefix}.6mA.bed.gz"
    File? ma_pileup_bed_index = "~{out_prefix}.6mA.bed.gz.tbi"
    File? summary = "~{out_prefix}.summary.json"
    String stat_hap1_cpg_count = read_string("~{out_prefix}.hap1.bed.count")
    String stat_hap2_cpg_count = read_string("~{out_prefix}.hap2.bed.count")
    String stat_combined_cpg_count = read_string("~{out_prefix}.combined.bed.count")
    Array[String] msg = read_lines("messages.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build2
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

task methbat_profile {
  meta {
    description: "Profile methylation regions with MethBat"
    outputs: {
      profile: {
        description: "MethBat 5mCpG profile"
      },
      stat_methbat_methylated_count: {
        description: "Number of profiled regions labeled as methylated"
      },
      stat_methbat_unmethylated_count: {
        description: "Number of profiled regions labeled as unmethylated"
      },
      stat_methbat_asm_count: {
        description: "Number of profiled regions labeled as having allele-specific methylation"
      }
    }
  }

  parameter_meta {
    cpg_pileup_bed: {
      description: "5mCpG pileup BED"
    }
    region_tsv: {
      description: "Regions TSV"
    }
    out_prefix: {
      description: "Output prefix"
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
    File cpg_pileup_bed
    File region_tsv
    String out_prefix
    Int threads = 2
    Int mem_gb = 4
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(cpg_pileup_bed, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    ln --symbolic --verbose "~{cpg_pileup_bed}" .

    methbat --version

    methbat profile \
      --input-pileup "~{basename(cpg_pileup_bed)}" \
      --input-regions "~{region_tsv}" \
      --output-region-profile "~{out_prefix}.methbat.profile.tsv"

    # count for three most interesting methylation summary_labels
    awk '$5=="Methylated" {print}' "~{out_prefix}.methbat.profile.tsv" \
    | wc -l > methylated_count.txt
    awk '$5=="Unmethylated" {print}' "~{out_prefix}.methbat.profile.tsv" \
    | wc -l > unmethylated_count.txt
    awk '$5=="AlleleSpecificMethylation" {print}' "~{out_prefix}.methbat.profile.tsv" \
    | wc -l > asm_count.txt
  >>>

  output {
    File profile = "~{out_prefix}.methbat.profile.tsv"
    String stat_methbat_methylated_count = read_string("methylated_count.txt")
    String stat_methbat_unmethylated_count = read_string("unmethylated_count.txt")
    String stat_methbat_asm_count = read_string("asm_count.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:281569947c0a6154f9a6bdaa1bb5f1e67dd160ccb49028d3218a49a37a0488fb"  # 1.1.0_build2
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

