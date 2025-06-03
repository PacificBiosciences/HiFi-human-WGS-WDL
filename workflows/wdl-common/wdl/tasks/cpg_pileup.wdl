version 1.0

import "../structs.wdl"

task cpg_pileup {
  meta {
    description: "Generate CpG methylation scores from aligned reads"
  }

  parameter_meta {
    haplotagged_bam: {
      name: "Aligned BAM"
    }
    haplotagged_bam_index: {
      name: "Aligned BAM index"
    }
    min_mapq: {
      name: "Minimum mapping quality"
    }
    min_coverage: {
      name: "Minimum coverage"
    }
    out_prefix: {
      name: "Output prefix"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    combined_bed: {
      name: "5mCpG BED for all alignments"
    }
    combined_bed_index: {
      name: "5mCpG BED index for all alignments"
    }
    hap1_bed: {
      name: "5mCpG BED for HP1 alignments"
    }
    hap1_bed_index: {
      name: "5mCpG BED index for HP1 alignments"
    }
    hap2_bed: {
      name: "5mCpG BED for HP2 alignments"
    }
    hap2_bed_index: {
      name: "5mCpG BED index for HP2 alignments"
    }
    combined_bw: {
      name: "5mCpG bigWig for all alignments"
    }
    hap1_bw: {
      name: "5mCpG bigWig for HP1 alignments"
    }
    hap2_bw: {
      name: "5mCpG bigWig for HP2 alignments"
    }
  }

  input {
    File haplotagged_bam
    File haplotagged_bam_index

    Int min_mapq     = 1
    Int min_coverage = 4

    String out_prefix

    File ref_fasta
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 12
  Int mem_gb    = 12
  Int disk_size = ceil((size(haplotagged_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -eu

    aligned_bam_to_cpg_scores --version

    aligned_bam_to_cpg_scores \
      --threads ~{threads} \
      --bam ~{haplotagged_bam} \
      --ref ~{ref_fasta} \
      --output-prefix ~{out_prefix} \
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
    File?  combined_bed            = "~{out_prefix}.cpg_pileup.combined.bed.gz"
    File?  combined_bed_index      = "~{out_prefix}.cpg_pileup.combined.bed.gz.tbi"
    File?  hap1_bed                = "~{out_prefix}.cpg_pileup.hap1.bed.gz"
    File?  hap1_bed_index          = "~{out_prefix}.cpg_pileup.hap1.bed.gz.tbi"
    File?  hap2_bed                = "~{out_prefix}.cpg_pileup.hap2.bed.gz"
    File?  hap2_bed_index          = "~{out_prefix}.cpg_pileup.hap2.bed.gz.tbi"
    File?  combined_bw             = "~{out_prefix}.cpg_pileup.combined.bw"
    File?  hap1_bw                 = "~{out_prefix}.cpg_pileup.hap1.bw"
    File?  hap2_bw                 = "~{out_prefix}.cpg_pileup.hap2.bw"
    String stat_hap1_cpg_count     = read_string("~{out_prefix}.hap1.bed.count")
    String stat_hap2_cpg_count     = read_string("~{out_prefix}.hap2.bed.count")
    String stat_combined_cpg_count = read_string("~{out_prefix}.combined.bed.count")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb-cpg-tools@sha256:afd5468a423fe089f1437d525fdc19c704296f723958739a6fe226caa01fba1c"
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
