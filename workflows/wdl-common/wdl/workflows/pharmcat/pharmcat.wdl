version 1.0

import "../../structs.wdl"

workflow pharmcat {
	input {
		IndexData haplotagged_bam
		IndexData phased_vcf

		IndexData reference

		String pangu_mode = "capture"

		IndexData pharmcat_positions
		Int pharmcat_min_coverage

		RuntimeAttributes default_runtime_attributes
	}

	call pangu_cyp2d6 {
		input:
			haplotagged_bam = haplotagged_bam.data,
			haplotagged_bam_index = haplotagged_bam.data_index,
			mode = pangu_mode,
			runtime_attributes = default_runtime_attributes
	}

	call pharmcat_preprocess {
		input:
			vcf = phased_vcf.data,
			vcf_index = phased_vcf.data_index,
			reference = reference.data,
			reference_index = reference.data_index,
			pharmcat_positions = pharmcat_positions.data,
			pharmcat_positions_index = pharmcat_positions.data_index,
			runtime_attributes = default_runtime_attributes
	}

	call filter_preprocessed_vcf {
		input:
			preprocessed_vcf = pharmcat_preprocess.preprocessed_vcf,
			haplotagged_bam = haplotagged_bam.data,
			haplotagged_bam_index = haplotagged_bam.data_index,
			reference_index = reference.data_index,
			min_coverage = pharmcat_min_coverage,
			runtime_attributes = default_runtime_attributes
	}

	call run_pharmcat {
		input:
			preprocessed_filtered_vcf = filter_preprocessed_vcf.filtered_vcf,
			pangu_tsv = pangu_cyp2d6.fixed_pangu_tsv,
			runtime_attributes = default_runtime_attributes
	}

	output {
		File pangu_json = pangu_cyp2d6.pangu_json
		File pangu_tsv = pangu_cyp2d6.pangu_tsv
		File fixed_pangu_tsv = pangu_cyp2d6.fixed_pangu_tsv

		File? pharmcat_missing_pgx_vcf = pharmcat_preprocess.missing_pgx_vcf
		File pharmcat_preprocessed_filtered_vcf = filter_preprocessed_vcf.filtered_vcf

		File pharmcat_match_json = run_pharmcat.pharmcat_match_json
		File pharmcat_phenotype_json = run_pharmcat.pharmcat_phenotype_json
		File pharmcat_report_html = run_pharmcat.pharmcat_report_html
		File pharmcat_report_json = run_pharmcat.pharmcat_report_json
	}

	parameter_meta {
		haplotagged_bam: {help: "Haplotagged BAM and index"}
		phased_vcf: {help: "Phased VCF and index"}
		reference: {help: "Reference genome fasta and index"}
		pangu_mode: {help: "Input type for pangu [\"wgs\", \"amplicon\", \"capture\" , \"consensus\"] (default = \"capture\")"}
		pharmcat_positions: {help: "VCF file and index specifying Pharmcat positions"}
		pharmcat_min_coverage: {help: "Minimum coverage cutoff used to filter the preprocessed VCF passed to pharmcat"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

# Call CYP2D6 for sample
task pangu_cyp2d6 {
	input {
		File haplotagged_bam
		File haplotagged_bam_index

		String mode

		RuntimeAttributes runtime_attributes
	}

	String haplotagged_bam_basename = basename(haplotagged_bam, ".bam")
	Int disk_size = ceil(size(haplotagged_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		pangu \
			-m ~{mode} \
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
		docker: "~{runtime_attributes.container_registry}/pangu@sha256:4d7d121280f6c30a6281ca7668530a5b129ee1f12cf561545caadeb2df1b83a8"
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

		File haplotagged_bam
		File haplotagged_bam_index

		File reference_index

		Int min_coverage

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(preprocessed_vcf, ".vcf.bgz")
	Int disk_size = ceil((size(preprocessed_vcf, "GB") + size(haplotagged_bam, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		bedtools coverage \
			-sorted \
			-g ~{reference_index} \
			-f 1 \
			-header \
			-mean \
			-a ~{preprocessed_vcf} \
			-b ~{haplotagged_bam} \
		| ( sed  -u '/^#CHROM/q' ; awk '$11 >= ~{min_coverage}' | cut -f1-10 ) \
		> ~{vcf_basename}.filtered.vcf
	>>>

	output {
		File filtered_vcf = "~{vcf_basename}.filtered.vcf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:cbe496e16773d4ad6f2eec4bd1b76ff142795d160f9dd418318f7162dcdaa685"
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

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(preprocessed_filtered_vcf, ".vcf")
	Int disk_size = ceil(size(preprocessed_filtered_vcf, "GB") * 2 + 20)

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
