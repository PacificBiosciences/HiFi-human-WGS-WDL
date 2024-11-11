version 1.0

import "../structs.wdl"

task hificnv {
  meta {
    description: "Call copy number reads from HiFi reads based on depth with HiFiCNV"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    vcf: {
      name: "Small variant VCF"
    }
    vcf_index: {
      name: "Small variant VCF index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    ref_name: {
      name: "Reference name"
    }
    exclude_bed: {
      name: "Regions to exclude from CNV calls"
    }
    exclude_bed_index: {
      name: "Regions to exclude from CNV calls (index)"
    }
    expected_male_bed: {
      name: "Expected CN BED for sample with XY karyotype"
    }
    expected_female_bed: {
      name: "Expected CN BED for sample with XX karyotype"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    cnv_vcf: {
      name: "HiFiCNV VCF"
    }
    cnv_vcf_index: {
      name: "HiFiCNV VCF index"
    }
    copynum_bedgraph: {
      name: "HiFiCNV CN bedgraph"
    }
    depth_bw: {
      name: "Depth bigWig"
    }
    maf_bw: {
      name: "MAF bigWig"
    }
    stat_DUP_count: {
      name: "Number of passing duplications"
    }
    stat_DUP_sum: {
      name: "Total length of passing duplications (Mbp)"
    }
    stat_DEL_count: {
      name: "Number of passing deletions"
    }
    stat_DEL_sum: {
      name: "Total length of passing deletions (Mbp)"
    }
  }

  input {
    String sample_id
    String? sex

    File aligned_bam
    File aligned_bam_index

    File vcf
    File vcf_index

    File ref_fasta
    File ref_index
    String ref_name

    File exclude_bed
    File exclude_bed_index

    File expected_male_bed
    File expected_female_bed

    RuntimeAttributes runtime_attributes
  }

  File expected_bed = if select_first([sex, "FEMALE"]) == "MALE" then expected_male_bed else expected_female_bed

  Int threads   = 8
  Int mem_gb    = threads * 2
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) + 20)

  command <<<
    set -euo pipefail

    echo ~{if defined(sex) then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for HiFiCNV."}

    hificnv --version
    bcftools --version

    hificnv \
      --threads ~{threads} \
      --bam ~{aligned_bam} \
      --ref ~{ref_fasta} \
      --maf ~{vcf} \
      --exclude ~{exclude_bed} \
      --expected-cn ~{expected_bed} \
      --output-prefix hificnv

    mv hificnv.~{sample_id}.copynum.bedgraph ~{sample_id}.~{ref_name}.hificnv.copynum.bedgraph
    mv hificnv.~{sample_id}.depth.bw ~{sample_id}.~{ref_name}.hificnv.depth.bw
    mv hificnv.~{sample_id}.maf.bw ~{sample_id}.~{ref_name}.hificnv.maf.bw
    mv hificnv.~{sample_id}.vcf.gz ~{sample_id}.~{ref_name}.hificnv.vcf.gz
    bcftools index --tbi ~{sample_id}.~{ref_name}.hificnv.vcf.gz

    # count the number of CNVs and sum the lengths of CNVs
    cat << EOF > cnv_stats.py
    import sys, pandas as pd
    df = pd.read_csv(sys.stdin, header=None)
    print(f'{df[0].count()}\t{df[0].sum()}')
    EOF

    bcftools query \
      --include 'FILTER="PASS" & SVTYPE="DUP"' \
      --format '%INFO/SVLEN\n' \
      ~{sample_id}.~{ref_name}.hificnv.vcf.gz \
      > DUP_lengths.txt
    if [ -s DUP_lengths.txt ]; then
      python3 ./cnv_stats.py < DUP_lengths.txt > stat_DUP.txt
      cut -f1 stat_DUP.txt > stat_DUP_count.txt
      cut -f2 stat_DUP.txt > stat_DUP_sum.txt
    else
      echo "0" > stat_DUP_count.txt
      echo "0" > stat_DUP_sum.txt
    fi

    bcftools query \
      --include 'FILTER="PASS" & SVTYPE="DEL"' \
      --format '%INFO/SVLEN\n' \
      ~{sample_id}.~{ref_name}.hificnv.vcf.gz \
      > DEL_lengths.txt
    if [ -s DEL_lengths.txt ]; then
      python3 ./cnv_stats.py < DEL_lengths.txt > stat_DEL.txt
      cut -f1 stat_DEL.txt > stat_DEL_count.txt
      cut -f2 stat_DEL.txt > stat_DEL_sum.txt
    else
      echo "0" > stat_DEL_count.txt
      echo "0" > stat_DEL_sum.txt
    fi
  >>>

  output {
    File   cnv_vcf          = "~{sample_id}.~{ref_name}.hificnv.vcf.gz"
    File   cnv_vcf_index    = "~{sample_id}.~{ref_name}.hificnv.vcf.gz.tbi"
    File   copynum_bedgraph = "~{sample_id}.~{ref_name}.hificnv.copynum.bedgraph"
    File   depth_bw         = "~{sample_id}.~{ref_name}.hificnv.depth.bw"
    File   maf_bw           = "~{sample_id}.~{ref_name}.hificnv.maf.bw"
    String stat_DUP_count   = read_string("stat_DUP_count.txt")
    String stat_DUP_sum     = read_string("stat_DUP_sum.txt")
    String stat_DEL_count   = read_string("stat_DEL_count.txt")
    String stat_DEL_sum     = read_string("stat_DEL_sum.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hificnv@sha256:c4764a70c8c2028edb1cdb4352997269947c5076ddd1aeaeef6c5076c630304d"
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
