version 1.0

import "../../structs.wdl"

workflow pharmcat {
  meta {
    description: "Run PharmCAT for a sample"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    phased_vcf: {
      name: "Phased small variant VCF"
    }
    phased_vcf_index: {
      name: "Phased small variant VCF index"
    }
    input_tsvs: {
      name: "StarPhase and HiFiHLA TSVs"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    pharmcat_version: {
      name: "PharmCAT Version"
    }
    pharmcat_positions: {
      name: "PharmCAT positions VCF"
    }
    pharmcat_positions_index: {
      name: "PharmCAT positions VCF index"
    }
    pharmcat_min_coverage: {
      name: "PharmCAT minimum coverage"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pharmcat_missing_pgx_vcf: {
      name: "PharmCAT missing PGx VCF"
    }
    pharmcat_preprocessed_filtered_vcf: {
      name: "PharmCAT preprocessed filtered VCF"
    }
    pharmcat_match_json: {
      name: "PharmCAT match JSON"
    }
    pharmcat_phenotype_json: {
      name: "PharmCAT phenotype JSON"
    }
    pharmcat_report_html: {
      name: "PharmCAT report HTML"
    }
    pharmcat_report_json: {
      name: "PharmCAT report JSON"
    }
  }

  input {
    String sample_id
    File haplotagged_bam
    File haplotagged_bam_index
    File phased_vcf
    File phased_vcf_index
    Array[File] input_tsvs

    File ref_fasta
    File ref_index

    String pharmcat_version
    File pharmcat_positions
    File pharmcat_positions_index
    Int pharmcat_min_coverage

    RuntimeAttributes default_runtime_attributes
  }

  String pharmcat_docker = (if (default_runtime_attributes.backend == "AWS-HealthOmics") then default_runtime_attributes.container_registry else "pgkb") + "/pharmcat:~{pharmcat_version}"

  call pharmcat_preprocess {
    input:
      vcf                      = phased_vcf,
      vcf_index                = phased_vcf_index,
      ref_fasta                = ref_fasta,
      ref_index                = ref_index,
      pharmcat_positions       = pharmcat_positions,
      pharmcat_positions_index = pharmcat_positions_index,
      pharmcat_docker          = pharmcat_docker,
      runtime_attributes       = default_runtime_attributes
  }

  call filter_preprocessed_vcf {
    input:
      preprocessed_vcf      = pharmcat_preprocess.preprocessed_vcf,
      haplotagged_bam       = haplotagged_bam,
      haplotagged_bam_index = haplotagged_bam_index,
      ref_index             = ref_index,
      min_coverage          = pharmcat_min_coverage,
      runtime_attributes    = default_runtime_attributes
  }

  if (defined(filter_preprocessed_vcf.filtered_vcf)) {
    call run_pharmcat {
      input:
        preprocessed_filtered_vcf = select_first([filter_preprocessed_vcf.filtered_vcf]),
        input_tsvs                = input_tsvs,
        pharmcat_docker           = pharmcat_docker,
        out_prefix                = "~{sample_id}.pharmcat",
        runtime_attributes        = default_runtime_attributes
    }
  }

  output {
    File? pharmcat_match_json     = run_pharmcat.pharmcat_match_json
    File? pharmcat_phenotype_json = run_pharmcat.pharmcat_phenotype_json
    File? pharmcat_report_html    = run_pharmcat.pharmcat_report_html
    File? pharmcat_report_json    = run_pharmcat.pharmcat_report_json
  }
}

task pharmcat_preprocess {
  meta {
    description: "Preprocess phased VCF for PharmCAT"
  }

  parameter_meta {
    vcf: {
      name: "Phased small variant VCF"
    }
    vcf_index: {
      name: "Phased small variant VCF index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    pharmcat_positions: {
      name: "PharmCAT positions VCF"
    }
    pharmcat_positions_index: {
      name: "PharmCAT positions VCF index"
    }
    pharmcat_docker: {
      name: "PharmCAT Docker image"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    preprocessed_vcf: {
      name: "Preprocessed VCF"
    }
    missing_pgx_vcf: {
      name: "Missing PGx VCF"
    }
  }

  input {
    File vcf
    File vcf_index

    File ref_fasta
    File ref_index

    File pharmcat_positions
    File pharmcat_positions_index

    String pharmcat_docker

    RuntimeAttributes runtime_attributes
  }

  String out_prefix = basename(vcf, ".vcf.gz")

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil((size(vcf, "GB") + size(ref_fasta, "GB") + size(pharmcat_positions, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    /pharmcat/pharmcat_vcf_preprocessor.py \
      --missing-to-ref \
      -vcf ~{vcf} \
      -refFna ~{ref_fasta} \
      -refVcf ~{pharmcat_positions} \
      -o .
  >>>

  output {
    File preprocessed_vcf = "~{out_prefix}.preprocessed.vcf.bgz"
    File? missing_pgx_vcf = "~{out_prefix}.missing_pgx_var.vcf"
  }

  runtime {
    docker: pharmcat_docker
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

# Remove ref calls with low mean coverage for sample
task filter_preprocessed_vcf {
  meta {
    description: "Filter preprocessed VCF for sample"
  }

  parameter_meta {
    preprocessed_vcf: {
      name: "Preprocessed VCF"
    }
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    min_coverage: {
      name: "Minimum coverage"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    filtered_vcf: {
      name: "Filtered VCF"
    }
  }

  input {
    File preprocessed_vcf

    File haplotagged_bam
    File haplotagged_bam_index

    File ref_index

    Int min_coverage

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 4
  Int disk_size = ceil((size(preprocessed_vcf, "GB") + size(haplotagged_bam, "GB")) + 20)

  String preprocessed_vcf_basename = basename(preprocessed_vcf)
  String out_prefix                = basename(preprocessed_vcf, ".vcf.bgz")

  command <<<
    set -euo pipefail

    bcftools --version

    ln --symbolic --verbose ~{preprocessed_vcf} .
    bcftools index --tbi ~{preprocessed_vcf_basename}

    # get a list of the regions overlapping variants
    bcftools query \
      --format '%CHROM\t%POS0\t%END\n' \
      ~{preprocessed_vcf} \
    > targeted.bed

    mosdepth --version

    # get the mean coverage for each region
    mosdepth \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --by targeted.bed \
      --no-per-base \
      targeted \
      ~{haplotagged_bam}

    # filter the bed for regions >= min_coverage
    cat << EOF > filter.py
    import gzip
    with gzip.open('targeted.regions.bed.gz', 'rt') as f:
      for line in f.readlines():
        if float(line.split('\t')[3]) >= ~{min_coverage}:
          print(line.strip())
    EOF

    python3 ./filter.py > targeted_regions.sufficient_depth.bed

    # filter the vcf for regions >= min_coverage
    if [ -s targeted_regions.sufficient_depth.bed ]; then
      bcftools view \
        --regions-file targeted_regions.sufficient_depth.bed \
        --output-type v \
        --output ~{out_prefix}.filtered.vcf \
        ~{preprocessed_vcf_basename}
    fi
  >>>

  output {
    File? filtered_vcf = "~{out_prefix}.filtered.vcf"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:63f7a5d1a4a17b71e66d755d3301a951e50f6b63777d34dab3ee9e182fd7acb1"
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

task run_pharmcat {
  meta {
    description: "Run PharmCAT for sample"
  }

  parameter_meta {
    preprocessed_filtered_vcf: {
      name: "Preprocessed filtered VCF"
    }
    input_tsvs: {
      name: "Pangu, StarPhase, and/or HiFiHLA TSVs"
    }
    pharmcat_docker: {
      name: "PharmCAT Docker image"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pharmcat_match_json: {
      name: "PharmCAT match JSON"
    }
    pharmcat_phenotype_json: {
      name: "PharmCAT phenotype JSON"
    }
    pharmcat_report_html: {
      name: "PharmCAT report HTML"
    }
    pharmcat_report_json: {
      name: "PharmCAT report JSON"
    }
  }

  input {
    File preprocessed_filtered_vcf
    Array[File] input_tsvs

    String pharmcat_docker

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(preprocessed_filtered_vcf, "GB") * 2 + 20)

  String pharmcat_basename = basename(preprocessed_filtered_vcf, ".vcf")

  command <<<
    set -eu

    sort -k1,1 < ~{sep=" " input_tsvs} \
    | grep -Ev 'NO_READS|NO_MATCH' \
    > merged.tsv \
    || touch merged.tsv

    # Run pharmcat
    /pharmcat/pharmcat \
      -vcf ~{preprocessed_filtered_vcf} \
      -reporterJson \
      "$([ -s merged.tsv ] && echo '-po merged.tsv')" \
      -o .
    
    # rename output files
    mv "~{pharmcat_basename}.match.json" "~{out_prefix}.match.json"
    mv "~{pharmcat_basename}.phenotype.json" "~{out_prefix}.phenotype.json"
    mv "~{pharmcat_basename}.report.html" "~{out_prefix}.report.html"
    mv "~{pharmcat_basename}.report.json" "~{out_prefix}.report.json"
  >>>

  output {
    File pharmcat_match_json     = "~{out_prefix}.match.json"
    File pharmcat_phenotype_json = "~{out_prefix}.phenotype.json"
    File pharmcat_report_html    = "~{out_prefix}.report.html"
    File pharmcat_report_json    = "~{out_prefix}.report.json"
  }

  runtime {
    docker: pharmcat_docker
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
