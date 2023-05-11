version 1.0

# Annotate small and structural variant VCFs using slivar. Outputs annotated VCFs and TSVs.
# This workflow is run on a phased single-sample VCF if there is only a single individual in the cohort, otherwise it is run on the joint-called phased VCF.

import "../humanwgs_structs.wdl"
import "../wdl-common/wdl/tasks/zip_index_vcf.wdl" as ZipIndexVcf

workflow tertiary_analysis {
	input {
		Cohort cohort
		IndexData small_variant_vcf
		IndexData sv_vcf

		ReferenceData reference

		SlivarData slivar_data

		RuntimeAttributes default_runtime_attributes
	}

	call write_cohort_yaml {
		input:
			cohort_id = cohort.cohort_id,
			cohort_json = write_json(cohort),
			runtime_attributes = default_runtime_attributes
	}

	call write_ped {
		input:
			cohort_id = cohort.cohort_id,
			cohort_yaml = write_cohort_yaml.cohort_yaml,
			runtime_attributes = default_runtime_attributes
	}

	call calculate_phrank {
		input:
			cohort_id = cohort.cohort_id,
			cohort_yaml = write_cohort_yaml.cohort_yaml,
			hpo_terms = slivar_data.hpo_terms,
			hpo_dag = slivar_data.hpo_dag,
			hpo_annotations = slivar_data.hpo_annotations,
			ensembl_to_hgnc = slivar_data.ensembl_to_hgnc,
			runtime_attributes = default_runtime_attributes
	}

	call bcftools_norm {
		input:
			vcf = small_variant_vcf.data,
			vcf_index = small_variant_vcf.data_index,
			reference = reference.fasta.data,
			runtime_attributes = default_runtime_attributes
	}

	call slivar_small_variant {
		input:
			bcf = bcftools_norm.normalized_bcf,
			bcf_index = bcftools_norm.normalized_bcf_index,
			pedigree = write_ped.pedigree,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			slivar_js = slivar_data.slivar_js,
			gnomad_af = reference.gnomad_af,
			hprc_af = reference.hprc_af,
			gff = reference.gff,
			runtime_attributes = default_runtime_attributes
	}

	call slivar_compound_hets {
		input:
			filtered_vcf = slivar_small_variant.filtered_vcf,
			filtered_vcf_index = slivar_small_variant.filtered_vcf_index,
			pedigree = write_ped.pedigree,
			runtime_attributes = default_runtime_attributes
	}

	call slivar_tsv {
		input:
			filtered_vcf = slivar_small_variant.filtered_vcf,
			compound_het_vcf = slivar_compound_hets.compound_het_vcf,
			pedigree = write_ped.pedigree,
			lof_lookup = slivar_data.lof_lookup,
			clinvar_lookup = slivar_data.clinvar_lookup,
			phrank_lookup = calculate_phrank.phrank_lookup,
			runtime_attributes = default_runtime_attributes
	}

	scatter (vcf_object in reference.population_vcfs) {
		File population_vcf = vcf_object.data
		File population_vcf_index = vcf_object.data_index
	}

	call svpack_filter_annotated {
		input:
			sv_vcf = sv_vcf.data,
			population_vcfs = population_vcf,
			population_vcf_indices = population_vcf_index,
			gff = reference.gff,
			runtime_attributes = default_runtime_attributes
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = svpack_filter_annotated.svpack_vcf,
			runtime_attributes = default_runtime_attributes
	}

	call slivar_svpack_tsv {
		input:
			filtered_vcf = zip_index_vcf.zipped_vcf,
			pedigree = write_ped.pedigree,
			lof_lookup = slivar_data.lof_lookup,
			clinvar_lookup = slivar_data.clinvar_lookup,
			phrank_lookup = calculate_phrank.phrank_lookup,
			runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData filtered_small_variant_vcf = {"data": slivar_small_variant.filtered_vcf, "data_index": slivar_small_variant.filtered_vcf_index}
		IndexData compound_het_small_variant_vcf = {"data": slivar_compound_hets.compound_het_vcf, "data_index": slivar_compound_hets.compound_het_vcf_index}
		File filtered_small_variant_tsv = slivar_tsv.filtered_tsv
		File compound_het_small_variant_tsv = slivar_tsv.compound_het_tsv
		IndexData filtered_svpack_vcf = {"data": zip_index_vcf.zipped_vcf, "data_index": zip_index_vcf.zipped_vcf_index}
		File filtered_svpack_tsv = slivar_svpack_tsv.svpack_tsv
	}

	parameter_meta {
		cohort: {help: "Sample information for the cohort"}
		small_variant_vcf: {help: "Small variant VCF to annotate using slivar"}
		sv_vcf: {help: "Structural variant VCF to annotate using slivar"}
		reference: {help: "Reference genome data"}
		slivar_data: {help: "Data files used for annotation with slivar"}
		default_runtime_attributes: {help: "Default RuntimeAttributes; spot if preemptible was set to true, otherwise on_demand"}
	}
}

task write_cohort_yaml {
	input {
		String cohort_id
		File cohort_json

		RuntimeAttributes runtime_attributes
	}

	command <<<
		set -euo pipefail

		parse_cohort.py \
			--cohort_json ~{cohort_json} \
			--write_cohort_yaml ~{cohort_id}.yml
	>>>

	output {
		File cohort_yaml = "~{cohort_id}.yml"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/parse-cohort@sha256:94444e7e3fd151936c9bbcb8a64b6a5e7d8c59de53b256a83f15c4ea203977b4"
		cpu: 2
		memory: "4 GB"
		disk: "20 GB"
		disks: "local-disk " + "20" + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task write_ped {
	input {
		String cohort_id
		File cohort_yaml

		RuntimeAttributes runtime_attributes
	}

	command <<<
		set -euo pipefail

		yaml2ped.py \
			~{cohort_yaml} \
			~{cohort_id} \
			~{cohort_id}.ped
	>>>

	output {
		File pedigree = "~{cohort_id}.ped"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pyyaml@sha256:7ce6b587e6a75df225eda677f47aba060ae98da10ff9a1fbe345a43e6ac3147b"
		cpu: 2
		memory: "4 GB"
		disk: "20 GB"
		disks: "local-disk " + "20" + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task calculate_phrank {
	input {
		String cohort_id
		File cohort_yaml

		File hpo_terms
		File hpo_dag
		File hpo_annotations
		File ensembl_to_hgnc

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(hpo_terms, "GB") + size(hpo_dag, "GB") + size(hpo_annotations, "GB") + size(ensembl_to_hgnc, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/calculate_phrank.py \
			~{hpo_terms} \
			~{hpo_dag} \
			~{hpo_annotations} \
			~{ensembl_to_hgnc} \
			~{cohort_yaml} \
			~{cohort_id} \
			~{cohort_id}_phrank.tsv
	>>>

	output {
		File phrank_lookup = "~{cohort_id}_phrank.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pyyaml@sha256:7ce6b587e6a75df225eda677f47aba060ae98da10ff9a1fbe345a43e6ac3147b"
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

task bcftools_norm {
	input {
		File vcf
		File vcf_index

		File reference

		RuntimeAttributes runtime_attributes
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int threads = 2
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools --version

		bcftools norm \
			--threads ~{threads - 1} \
			--multiallelics \
			- \
			--output-type b \
			--fasta-ref ~{reference} \
			~{vcf} \
		| bcftools sort \
			--output-type b \
			--output ~{vcf_basename}.norm.bcf

		bcftools index --threads ~{threads - 1} ~{vcf_basename}.norm.bcf
	>>>

	output {
		File normalized_bcf = "~{vcf_basename}.norm.bcf"
		File normalized_bcf_index = "~{vcf_basename}.norm.bcf.csi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
		cpu: threads
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

task slivar_small_variant {
	input {
		File bcf
		File bcf_index
		File pedigree

		File reference
		File reference_index

		File slivar_js
		File gnomad_af
		File hprc_af
		File gff

		RuntimeAttributes runtime_attributes
	}

	Float max_gnomad_af = 0.01
	Float max_hprc_af = 0.01
	Int max_gnomad_nhomalt = 4
	Int max_hprc_nhomalt = 4
	Int max_gnomad_ac = 4
	Int max_hprc_ac = 4
	Int min_gq = 5

	Array[String] info_expr = [
		'variant.FILTER=="PASS"',
		'INFO.gnomad_af <= ~{max_gnomad_af}',
		'INFO.hprc_af <= ~{max_hprc_af}',
		'INFO.gnomad_nhomalt <= ~{max_gnomad_nhomalt}',
		'INFO.hprc_nhomalt <= ~{max_hprc_nhomalt}'
	]
	Array[String] family_recessive_expr = [
		'recessive:fam.every(segregating_recessive)'
	]
	Array[String] family_x_recessive_expr = [
		'x_recessive:(variant.CHROM == "chrX")',
		'fam.every(segregating_recessive_x)'
	]
	Array[String] family_dominant_expr = [
		'dominant:fam.every(segregating_dominant)',
		'INFO.gnomad_ac <= ~{max_gnomad_ac}',
		'INFO.hprc_ac <= ~{max_hprc_ac}'
	]
	Array[String] sample_expr = [
		'comphet_side:sample.het',
		'sample.GQ > ~{min_gq}'
	]

	String bcf_basename = basename(bcf, ".bcf")
	Int threads = 8
	Int disk_size = ceil((size(bcf, "GB") + size(reference, "GB") + size(gnomad_af, "GB") + size(hprc_af, "GB") + size(gff, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		# slivar has no version option
		slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+' 

		bcftools --version
		
		tabix --version

		pslivar \
			--processes ~{threads} \
			--fasta ~{reference} \
			--pass-only \
			--js ~{slivar_js} \
			--info '~{sep=" && " info_expr}' \
			--family-expr '~{sep=" && " family_recessive_expr}' \
			--family-expr '~{sep=" && " family_x_recessive_expr}' \
			--family-expr '~{sep=" && " family_dominant_expr}' \
			--sample-expr '~{sep=" && " sample_expr}' \
			--gnotate ~{gnomad_af} \
			--gnotate ~{hprc_af} \
			--vcf ~{bcf} \
			--ped ~{pedigree} \
		| bcftools csq \
			--local-csq \
			--samples - \
			--ncsq 40 \
			--gff-annot ~{gff} \
			--fasta-ref ~{reference} \
			- \
			--output-type z \
			--output ~{bcf_basename}.slivar.vcf.gz

		tabix ~{bcf_basename}.slivar.vcf.gz
	>>>

	output {
		File filtered_vcf = "~{bcf_basename}.slivar.vcf.gz"
		File filtered_vcf_index = "~{bcf_basename}.slivar.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/slivar@sha256:1148503061ca622fcdda1b280a3b13ccc502728140775338689f5cfdf89c4556"
		cpu: threads
		memory: "16 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task slivar_compound_hets {
	input {
		File filtered_vcf
		File filtered_vcf_index
		File pedigree

		RuntimeAttributes runtime_attributes
	}

	Array[String] skip_list = [
		'non_coding_transcript',
		'intron',
		'non_coding',
		'upstream_gene',
		'downstream_gene',
		'non_coding_transcript_exon',
		'NMD_transcript',
		'5_prime_UTR',
		'3_prime_UTR'
	]

	String vcf_basename = basename(filtered_vcf, ".vcf.gz")
	Int disk_size = ceil(size(filtered_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		# slivar has no version option
		slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+'

		bgzip --version

		tabix --version

		slivar \
			compound-hets \
			--skip ~{sep=',' skip_list} \
			--vcf ~{filtered_vcf} \
			--sample-field comphet_side \
			--ped ~{pedigree} \
			--allow-non-trios \
		| add_comphet_phase.py \
		> ~{vcf_basename}.compound_hets.vcf

		bgzip ~{vcf_basename}.compound_hets.vcf
		tabix ~{vcf_basename}.compound_hets.vcf.gz
	>>>

	output {
		File compound_het_vcf = "~{vcf_basename}.compound_hets.vcf.gz"
		File compound_het_vcf_index = "~{vcf_basename}.compound_hets.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/slivar@sha256:1148503061ca622fcdda1b280a3b13ccc502728140775338689f5cfdf89c4556"
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

task slivar_tsv {
	input {
		File filtered_vcf
		File compound_het_vcf
		File pedigree

		File lof_lookup
		File clinvar_lookup
		File phrank_lookup

		RuntimeAttributes runtime_attributes
	}

	Array[String] info_fields = [
		'gnomad_af',
		'hprc_af',
		'gnomad_nhomalt',
		'hprc_nhomalt',
		'gnomad_ac',
		'hprc_ac'
	]

	String filtered_vcf_basename = basename(filtered_vcf, ".vcf.gz")
	String compound_het_vcf_basename = basename(compound_het_vcf, ".vcf.gz")
	Int disk_size = ceil((size(filtered_vcf, "GB") + size(compound_het_vcf, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		# slivar has no version option
		slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+'

		slivar tsv \
			--info-field ~{sep=' --info-field ' info_fields} \
			--sample-field dominant \
			--sample-field recessive \
			--sample-field x_recessive \
			--csq-field BCSQ \
			--gene-description ~{lof_lookup} \
			--gene-description ~{clinvar_lookup} \
			--gene-description ~{phrank_lookup} \
			--ped ~{pedigree} \
			--out /dev/stdout \
			~{filtered_vcf} \
		| sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
		> ~{filtered_vcf_basename}.tsv

		slivar tsv \
			--info-field ~{sep=' --info-field ' info_fields} \
			--sample-field slivar_comphet \
			--info-field slivar_comphet \
			--csq-field BCSQ \
			--gene-description ~{lof_lookup} \
			--gene-description ~{clinvar_lookup} \
			--gene-description ~{phrank_lookup} \
			--ped ~{pedigree} \
			--out /dev/stdout \
			~{compound_het_vcf} \
		| sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
		> ~{compound_het_vcf_basename}.tsv
	>>>

	output {
		File filtered_tsv = "~{filtered_vcf_basename}.tsv"
		File compound_het_tsv = "~{compound_het_vcf_basename}.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/slivar@sha256:1148503061ca622fcdda1b280a3b13ccc502728140775338689f5cfdf89c4556"
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

task svpack_filter_annotated {
	input {
		File sv_vcf

		Array[File] population_vcfs
		Array[File] population_vcf_indices

		File gff

		RuntimeAttributes runtime_attributes
	}

	String sv_vcf_basename = basename(sv_vcf, ".vcf.gz")
	Int disk_size = ceil(size(sv_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		echo "svpack version:"
		cat /opt/svpack/.git/HEAD

		svpack \
			filter \
			--pass-only \
			--min-svlen 50 \
			~{sv_vcf} \
		~{sep=' ' prefix('| svpack match -v - ', population_vcfs)} \
		| svpack \
			consequence \
			- \
			~{gff} \
		| svpack \
			tagzygosity \
			- \
		> ~{sv_vcf_basename}.svpack.vcf
	>>>

	output {
		File svpack_vcf = "~{sv_vcf_basename}.svpack.vcf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/svpack@sha256:cd357d4a032e566946d651b8870658cdeed7ee98f46b6e202a8b6e3edf553507"
		cpu: 2
		memory: "16 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}

task slivar_svpack_tsv {
	input {
		File filtered_vcf

		File pedigree
		File lof_lookup
		File clinvar_lookup
		File phrank_lookup

		RuntimeAttributes runtime_attributes
	}

	Array[String] info_fields = [
		'SVTYPE',
		'SVLEN',
		'SVANN',
		'CIPOS',
		'MATEID',
		'END'
	]

	String filtered_vcf_basename = basename(filtered_vcf, ".vcf.gz")
	Int disk_size = ceil((size(filtered_vcf, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		# slivar has no version option
		slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+'

		slivar tsv \
			--info-field ~{sep=' --info-field ' info_fields} \
			--sample-field hetalt \
			--sample-field homalt \
			--csq-field BCSQ \
			--gene-description ~{lof_lookup} \
			--gene-description ~{clinvar_lookup} \
			--gene-description ~{phrank_lookup} \
			--ped ~{pedigree} \
			--out /dev/stdout \
			~{filtered_vcf} \
		| sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
		> ~{filtered_vcf_basename}.tsv
	>>>

	output {
		File svpack_tsv = "~{filtered_vcf_basename}.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/slivar@sha256:1148503061ca622fcdda1b280a3b13ccc502728140775338689f5cfdf89c4556"
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
