version 1.0

import "../structs.wdl"

task cpg_pileup {
  meta {
    description: "Generate CpG methylation scores from aligned reads"
    outputs: {
      combined_bed: {
        description: "5mCpG combined BED"
      },
      combined_bed_index: {
        description: "Index for 5mCpG combined BED"
      },
      hap1_bed: {
        description: "5mCpG haplotype 1 BED"
      },
      hap1_bed_index: {
        description: "Index for 5mCpG haplotype 1 BED"
      },
      hap2_bed: {
        description: "5mCpG haplotype 2 BED"
      },
      hap2_bed_index: {
        description: "Index for 5mCpG haplotype 2 BED"
      },
      combined_bw: {
        description: "5mCpG combined BigWig"
      },
      hap1_bw: {
        description: "5mCpG haplotype 1 BigWig"
      },
      hap2_bw: {
        description: "5mCpG haplotype 2 BigWig"
      },
      stat_hap1_cpg_count: {
        description: "Number of scored reference 5mCpGs in haplotype 1"
      },
      stat_hap2_cpg_count: {
        description: "Number of scored reference 5mCpGs in haplotype 2"
      },
      stat_combined_cpg_count: {
        description: "Number of scored reference 5mCpGs combined"
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
      description: "Minimum mapping quality"
    }
    min_coverage: {
      description: "Minimum coverage"
    }
    out_prefix: {
      description: "Output prefix"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
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
    String out_prefix
    File ref_fasta
    File ref_index
    RuntimeAttributes runtime_attributes
  }

  Int threads = 16
  Int mem_gb = 32
  Int disk_size = ceil((size(haplotagged_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -eu

    ln --symbolic --verbose "~{haplotagged_bam}" .
    ln --symbolic --verbose "~{haplotagged_bam_index}" .
    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    aligned_bam_to_cpg_scores --version

    aligned_bam_to_cpg_scores \
      --threads ~{threads} \
      --bam "~{basename(haplotagged_bam)}" \
      --ref "~{basename(ref_fasta)}" \
      --output-prefix "~{out_prefix}" \
      --min-mapq ~{min_mapq} \
      --min-coverage ~{min_coverage}

    # rename for clarity
    for infix in combined hap1 hap2; do
      echo "0" > "~{out_prefix}.${infix}.bed.count"
      if [ -f "~{out_prefix}.${infix}.bw" ]; then
        mv --verbose "~{out_prefix}.${infix}.bw" "~{out_prefix}.cpg_pileup.${infix}.bw"
      fi
      if [ -f "~{out_prefix}.${infix}.bed.gz" ]; then
        zgrep --count --invert-match "^#" "~{out_prefix}.${infix}.bed.gz" \
          > "~{out_prefix}.${infix}.bed.count" || true  # ignore files with no lines
        mv --verbose "~{out_prefix}.${infix}.bed.gz" "~{out_prefix}.cpg_pileup.${infix}.bed.gz"
        mv --verbose "~{out_prefix}.${infix}.bed.gz.tbi" "~{out_prefix}.cpg_pileup.${infix}.bed.gz.tbi"
      fi
    done
  >>>

  output {
    File? combined_bed = "~{out_prefix}.cpg_pileup.combined.bed.gz"
    File? combined_bed_index = "~{out_prefix}.cpg_pileup.combined.bed.gz.tbi"
    File? hap1_bed = "~{out_prefix}.cpg_pileup.hap1.bed.gz"
    File? hap1_bed_index = "~{out_prefix}.cpg_pileup.hap1.bed.gz.tbi"
    File? hap2_bed = "~{out_prefix}.cpg_pileup.hap2.bed.gz"
    File? hap2_bed_index = "~{out_prefix}.cpg_pileup.hap2.bed.gz.tbi"
    File? combined_bw = "~{out_prefix}.cpg_pileup.combined.bw"
    File? hap1_bw = "~{out_prefix}.cpg_pileup.hap1.bw"
    File? hap2_bw = "~{out_prefix}.cpg_pileup.hap2.bw"
    String stat_hap1_cpg_count = read_string("~{out_prefix}.hap1.bed.count")
    String stat_hap2_cpg_count = read_string("~{out_prefix}.hap2.bed.count")
    String stat_combined_cpg_count = read_string("~{out_prefix}.combined.bed.count")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb-cpg-tools@sha256:afd5468a423fe089f1437d525fdc19c704296f723958739a6fe226caa01fba1c"  # 3.0.0_build1
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
