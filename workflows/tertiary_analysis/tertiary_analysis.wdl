version 1.0

import "../common/structs.wdl"
import "../common/tasks/zip_index_vcf.wdl" as ZipIndexVcf

workflow tertiary_analysis {
	input {
		Cohort cohort
		IndexData small_variant_vcf
		IndexData sv_vcf

		ReferenceData reference

		SlivarData slivar_data

		String container_registry
		Boolean preemptible
	}

	call write_cohort_yaml {
		input:
			cohort_id = cohort.cohort_id,
			cohort_json = write_json(cohort),
			container_registry = container_registry,
			preemptible = preemptible
	}

	call write_ped {
		input:
			cohort_id = cohort.cohort_id,
			cohort_yaml = write_cohort_yaml.cohort_yaml,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call calculate_phrank {
		input:
			cohort_id = cohort.cohort_id,
			cohort_yaml = write_cohort_yaml.cohort_yaml,
			hpo_terms = slivar_data.hpo_terms,
			hpo_dag = slivar_data.hpo_dag,
			hpo_annotations = slivar_data.hpo_annotations,
			ensembl_to_hgnc = slivar_data.ensembl_to_hgnc,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call bcftools_norm {
		input:
			vcf = small_variant_vcf.data,
			vcf_index = small_variant_vcf.data_index,
			reference = reference.fasta.data,
			container_registry = container_registry,
			preemptible = preemptible
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
			container_registry = container_registry,
			preemptible = preemptible
	}

	call slivar_compound_hets {
		input:
			filtered_vcf = slivar_small_variant.filtered_vcf,
			filtered_vcf_index = slivar_small_variant.filtered_vcf_index,
			pedigree = write_ped.pedigree,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call slivar_tsv {
		input:
			filtered_vcf = slivar_small_variant.filtered_vcf,
			compound_het_vcf = slivar_compound_hets.compound_het_vcf,
			pedigree = write_ped.pedigree,
			lof_lookup = slivar_data.lof_lookup,
			clinvar_lookup = slivar_data.clinvar_lookup,
			phrank_lookup = calculate_phrank.phrank_lookup,
			container_registry = container_registry,
			preemptible = preemptible
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
			container_registry = container_registry,
			preemptible = preemptible
	}

	call ZipIndexVcf.zip_index_vcf {
		input:
			vcf = svpack_filter_annotated.svpack_vcf,
			container_registry = container_registry,
			preemptible = preemptible
	}

	call slivar_svpack_tsv {
		input:
			filtered_vcf = zip_index_vcf.zipped_vcf,
			pedigree = write_ped.pedigree,
			lof_lookup = slivar_data.lof_lookup,
			clinvar_lookup = slivar_data.clinvar_lookup,
			phrank_lookup = calculate_phrank.phrank_lookup,
			container_registry = container_registry,
			preemptible = preemptible
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
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task write_cohort_yaml {
	input {
		String cohort_id
		File cohort_json

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/parse_cohort:1.0.0"
		cpu: 4
		memory: "14 GB"
		disk: "20 GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task write_ped {
	input {
		String cohort_id
		File cohort_yaml

		String container_registry
		Boolean preemptible
	}

	command <<<
		set -euo pipefail

		python3 /opt/scripts/yaml2ped.py \
			~{cohort_yaml} \
			~{cohort_id} \
			~{cohort_id}.ped
	>>>

	output {
		File pedigree = "~{cohort_id}.ped"
	}

	runtime {
		docker: "~{container_registry}/pyyaml:b1a46c6"
		cpu: 4
		memory: "14 GB"
		disk: "20 GB"
		preemptible: preemptible
		maxRetries: 3
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

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/pyyaml:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task bcftools_norm {
	input {
		File vcf
		File vcf_index

		File reference

		String container_registry
		Boolean preemptible
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools norm \
			--multiallelics \
			- \
			--output-type b \
			--fasta-ref ~{reference} \
			~{vcf} \
		| bcftools sort \
			--output-type b \
			--output ~{vcf_basename}.norm.bcf

		bcftools index ~{vcf_basename}.norm.bcf
	>>>

	output {
		File normalized_bcf = "~{vcf_basename}.norm.bcf"
		File normalized_bcf_index = "~{vcf_basename}.norm.bcf.csi"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
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

		String container_registry
		Boolean preemptible
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
	Array[String] family_x_dominant_expr = [
		'x_dominant:(variant.CHROM == "chrX")',
		'fam.every(segregating_dominant_x)',
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

		pslivar \
			--processes ~{threads} \
			--fasta ~{reference} \
			--pass-only \
			--js ~{slivar_js} \
			--info '~{sep=" && " info_expr}' \
			--family-expr '~{sep=" && " family_recessive_expr}' \
			--family-expr '~{sep=" && " family_x_recessive_expr}' \
			--family-expr '~{sep=" && " family_dominant_expr}' \
			--family-expr '~{sep=" && " family_x_dominant_expr}' \
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
		docker: "~{container_registry}/slivar:b1a46c6"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task slivar_compound_hets {
	input {
		File filtered_vcf
		File filtered_vcf_index
		File pedigree

		String container_registry
		Boolean preemptible
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

		slivar \
			compound-hets \
			--skip ~{sep=',' skip_list} \
			--vcf ~{filtered_vcf} \
			--sample-field comphet_side \
			--ped ~{pedigree} \
			--allow-non-trios \
		| python3 /opt/scripts/add_comphet_phase.py \
		> ~{vcf_basename}.compound_hets.vcf

		bgzip ~{vcf_basename}.compound_hets.vcf
		tabix ~{vcf_basename}.compound_hets.vcf.gz
	>>>

	output {
		File compound_het_vcf = "~{vcf_basename}.compound_hets.vcf.gz"
		File compound_het_vcf_index = "~{vcf_basename}.compound_hets.vcf.gz.tbi"
	}

	runtime {
		docker: "~{container_registry}/slivar:b1a46c6"
		cpu: 2
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
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

		String container_registry
		Boolean preemptible
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

		slivar tsv \
			--info-field ~{sep=' --info-field ' info_fields} \
			--sample-field dominant \
			--sample-field x_dominant \
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
		docker: "~{container_registry}/slivar:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task svpack_filter_annotated {
	input {
		File sv_vcf

		Array[File] population_vcfs
		Array[File] population_vcf_indices

		File gff

		String container_registry
		Boolean preemptible
	}

	String sv_vcf_basename = basename(sv_vcf, ".vcf.gz")
	Int disk_size = ceil(size(sv_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python /opt/scripts/svpack/svpack \
			filter \
			--pass-only \
			--min-svlen 50 \
			~{sv_vcf} \
		~{sep=' ' prefix('| python /opt/scripts/svpack/svpack match -v - ', population_vcfs)} \
		| python /opt/scripts/svpack/svpack \
			consequence \
			- \
			~{gff} \
		| python /opt/scripts/svpack/svpack \
			tagzygosity \
			- \
		> ~{sv_vcf_basename}.svpack.vcf
	>>>

	output {
		File svpack_vcf = "~{sv_vcf_basename}.svpack.vcf"
	}

	runtime {
		docker: "~{container_registry}/svpack:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}

task slivar_svpack_tsv {
	input {
		File filtered_vcf

		File pedigree
		File lof_lookup
		File clinvar_lookup
		File phrank_lookup

		String container_registry
		Boolean preemptible
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
		docker: "~{container_registry}/slivar:b1a46c6"
		cpu: 1
		memory: "1 GB"
		disk: disk_size + " GB"
		preemptible: preemptible
		maxRetries: 3
	}
}
