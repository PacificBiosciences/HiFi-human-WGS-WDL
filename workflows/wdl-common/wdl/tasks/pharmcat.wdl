version 1.0

import "../structs.wdl"

workflow pharmcat {
    input {
        Array[Pair[String,Map[String,IndexData]]] sample_data

        IndexData reference
        File reference_chromosome_lengths

        IndexData pharmcat_positions
        Int pharmcat_min_coverage

        RuntimeAttributes default_runtime_attributes
    }

    scatter (sample in sample_data) {
        call pangu_cyp2d6 {
            input:
                haplotagged_bam = sample.right["haplotagged_bam"].data,
                haplotagged_bam_index = sample.right["haplotagged_bam"].data_index,
                runtime_attributes = default_runtime_attributes
        }

        call pharmcat_preprocess {
            input:
                vcf = sample.right["phased_vcf"].data,
                vcf_index = sample.right["phased_vcf"].data_index,
                reference = reference.data,
                reference_index = reference.data_index,
                pharmcat_positions = pharmcat_positions.data,
                pharmcat_positions_index = pharmcat_positions.data_index,
                runtime_attributes = default_runtime_attributes
        }

        call filter_preprocessed_vcf {
            input:
                preprocessed_vcf = pharmcat_preprocess.preprocessed_vcf,
                aligned_bam = sample.right["aligned_bam"].data,
                aligned_bam_index = sample.right["aligned_bam"].data_index,
                reference_chromosome_lengths = reference_chromosome_lengths,
                min_coverage = pharmcat_min_coverage,
                runtime_attributes = default_runtime_attributes
        }

        call run_pharmcat {
            input:
                preprocessed_filtered_vcf = filter_preprocessed_vcf.filtered_vcf,
                pangu_tsv = pangu_cyp2d6.fixed_pangu_tsv,
                reference = reference.data,
                reference_index = reference.data_index,
                runtime_attributes = default_runtime_attributes
        }
    }

    output {
        Array[File] pangu_jsons = pangu_cyp2d6.pangu_json
        Array[File] pangu_tsvs = pangu_cyp2d6.pangu_tsv
        Array[File] fixed_pangu_tsvs = pangu_cyp2d6.fixed_pangu_tsv

        Array[File?] pharmcat_missing_pgx_vcfs = pharmcat_preprocess.missing_pgx_vcf
        Array[File] pharmcat_preprocessed_filtered_vcfs = filter_preprocessed_vcf.filtered_vcf

        Array[File] pharmcat_match_jsons = run_pharmcat.pharmcat_match_json
        Array[File] pharmcat_phenotype_jsons = run_pharmcat.pharmcat_phenotype_json
        Array[File] pharmcat_report_htmls = run_pharmcat.pharmcat_report_html
        Array[File] pharmcat_report_jsons = run_pharmcat.pharmcat_report_json
    }

    parameter_meta {
        sample_data: {help: "Array of pairs mapping sample ID to aligned bam, haplotagged bam, gvcf, and phased VCF files for the sample"}
        reference: {help: "Reference genome data"}
        pharmcat_positions: {help: "VCF file and index specifying pharmact positions"}
        pharmcat_min_coverage: {help: "Minimum coverage cutoff used to filter the preprocessed VCF passed to pharmcat"}
        default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
    }
}

# Call CYP2D6 for sample
task pangu_cyp2d6 {
    input {
        File haplotagged_bam
        File haplotagged_bam_index

        RuntimeAttributes runtime_attributes
    }

    String haplotagged_bam_basename = basename(haplotagged_bam, ".bam")
    Int disk_size = ceil(size(haplotagged_bam, "GB") * 2 + 20)

    command <<<
        set -euo pipefail

        pangu \
            -m capture \
            -p ~{haplotagged_bam_basename}.pangu \
            ~{haplotagged_bam}

        # Fix the pangu output with missing calls for the sample
        awk \
            'BEGIN {{OFS="\t"}} !($2 ~ /\//) {{$2=$2"/[]"}} 1' \
            ~{haplotagged_bam_basename}.pangu_pharmcat.tsv \
        > ~{haplotagged_bam_basename}.pangu_pharmcat_fix.tsv
    >>>

    output {
        File pangu_json = "~{haplotagged_bam_basename}.pangu_report.json"
        File pangu_tsv = "~{haplotagged_bam_basename}.pangu_pharmcat.tsv"
        File fixed_pangu_tsv = "~{haplotagged_bam_basename}.pangu_pharmcat_fix.tsv"
    }

    runtime {
        docker: "~{runtime_attributes.container_registry}/pangu@sha256:477dfa87eb98f54708dad3b20cab24ea1a171886b0b2b9d436b3ffc4e899b908"
        cpu: 2
        memory: "12 GB"
        disk: disk_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: runtime_attributes.preemptible_tries
        maxRetries: runtime_attributes.max_retries
        awsBatchRetryAttempts: runtime_attributes.max_retries
        queueArn: runtime_attributes.queue_arn
        zones: runtime_attributes.zones
    }
}

# Preprocess phased VCF for sample
task pharmcat_preprocess {
    input {
        File vcf
        File vcf_index

        File reference
        File reference_index

        File pharmcat_positions
        File pharmcat_positions_index

        RuntimeAttributes runtime_attributes
    }

    String vcf_basename = basename(vcf, ".vcf.gz")
    Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(pharmcat_positions, "GB")) * 2 + 20)

    command <<<
        set -euo pipefail

        /pharmcat/pharmcat_vcf_preprocessor.py \
            --missing-to-ref \
            -vcf ~{vcf} \
            -refFna ~{reference} \
            -refVcf ~{pharmcat_positions} \
            -o .
    >>>

    output {
        File preprocessed_vcf = "~{vcf_basename}.preprocessed.vcf.bgz"
        File? missing_pgx_vcf = "~{vcf_basename}.missing_pgx_var.vcf"
    }

    runtime {
        docker: "pgkb/pharmcat:2.3.0"
        cpu: 2
        memory: "4 GB"
        disk: disk_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: runtime_attributes.preemptible_tries
        maxRetries: runtime_attributes.max_retries
        awsBatchRetryAttempts: runtime_attributes.max_retries
        queueArn: runtime_attributes.queue_arn
        zones: runtime_attributes.zones
    }
}

# Remove ref calls with low mean coverage for sample
task filter_preprocessed_vcf {
    input {
        File preprocessed_vcf

        File aligned_bam
        File aligned_bam_index

        File reference_chromosome_lengths

        Int min_coverage

        RuntimeAttributes runtime_attributes
    }

    String vcf_basename = basename(preprocessed_vcf, ".vcf.bgz")
    Int disk_size = ceil((size(preprocessed_vcf, "GB") + size(aligned_bam, "GB")) * 2 + 20)

    command <<<
        set -euo pipefail

        bedtools coverage \
            -sorted \
            -g ~{reference_chromosome_lengths} \
            -f 1 \
            -header \
            -mean \
            -a ~{preprocessed_vcf} \
            -b ~{aligned_bam} \
        | ( sed  -u '/^#CHROM/q' ; awk '$11 >= ~{min_coverage}' | cut -f1-10 ) \
        > ~{vcf_basename}.filtered.vcf
    >>>

    output {
        File filtered_vcf = "~{vcf_basename}.filtered.vcf"
    }

    runtime {
        docker: "~{runtime_attributes.container_registry}/samtools@sha256:a843074b9be9505e6e6e93385975f761617fcce4c486fcebf97ab65075ed6bd4"
        cpu: 2
        memory: "4 GB"
        disk: disk_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: runtime_attributes.preemptible_tries
        maxRetries: runtime_attributes.max_retries
        awsBatchRetryAttempts: runtime_attributes.max_retries
        queueArn: runtime_attributes.queue_arn
        zones: runtime_attributes.zones
    }
}

# Run pharmcat for sample
task run_pharmcat {
    input {
        File preprocessed_filtered_vcf

        File pangu_tsv

        File reference
        File reference_index

        RuntimeAttributes runtime_attributes
    }

    String vcf_basename = basename(preprocessed_filtered_vcf, ".vcf")
    Int disk_size = ceil((size(preprocessed_filtered_vcf, "GB") + size(reference, "GB")) * 2 + 20)

    command <<<
        set -euo pipefail

        # Run pharmcat
        /pharmcat/pharmcat \
            -vcf ~{preprocessed_filtered_vcf} \
            -reporterJson \
            -po ~{pangu_tsv} \
            -o .
    >>>

    output {
        File pharmcat_match_json = "~{vcf_basename}.match.json"
        File pharmcat_phenotype_json = "~{vcf_basename}.phenotype.json"
        File pharmcat_report_html = "~{vcf_basename}.report.html"
        File pharmcat_report_json = "~{vcf_basename}.report.json"
    }

    runtime {
        docker: "pgkb/pharmcat:2.3.0"
        cpu: 2
        memory: "4 GB"
        disk: disk_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: runtime_attributes.preemptible_tries
        maxRetries: runtime_attributes.max_retries
        awsBatchRetryAttempts: runtime_attributes.max_retries
        queueArn: runtime_attributes.queue_arn
        zones: runtime_attributes.zones
    }
}
