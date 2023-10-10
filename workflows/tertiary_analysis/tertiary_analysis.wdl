version 1.0

# Annotate small and structural variant VCFs using slivar. Outputs annotated VCFs and TSVs.
# This workflow is run on a phased single-sample VCF if there is only a single individual in the cohort, otherwise it is run on the joint-called phased VCF.

import "../humanwgs_structs.wdl"

workflow tertiary_analysis {
	input {
		Cohort cohort
		IndexData small_variant_vcf
		IndexData sv_vcf

		ReferenceData reference

		SlivarData slivar_data

		RuntimeAttributes default_runtime_attributes
	}

	call write_yaml_ped_phrank {
		input:
			cohort_id = cohort.cohort_id,
			cohort_json = write_json(cohort),
			hpo_terms = slivar_data.hpo_terms,
			hpo_dag = slivar_data.hpo_dag,
			hpo_annotations = slivar_data.hpo_annotations,
			ensembl_to_hgnc = slivar_data.ensembl_to_hgnc,
			runtime_attributes = default_runtime_attributes
	}

	call slivar_small_variant {
		input:
			vcf = small_variant_vcf.data,
			vcf_index = small_variant_vcf.data_index,
			pedigree = write_yaml_ped_phrank.pedigree,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			slivar_js = slivar_data.slivar_js,
			gnomad_af = select_first([reference.gnomad_af]),
			hprc_af = select_first([reference.hprc_af]),
			gff = select_first([reference.gff]),
			lof_lookup = slivar_data.lof_lookup,
			clinvar_lookup = slivar_data.clinvar_lookup,
			phrank_lookup = write_yaml_ped_phrank.phrank_lookup,
			runtime_attributes = default_runtime_attributes
	}

	scatter (vcf_object in select_first([reference.population_vcfs])) {
		File population_vcf = vcf_object.data
		File population_vcf_index = vcf_object.data_index
	}

	call svpack_filter_annotated {
		input:
			sv_vcf = sv_vcf.data,
			population_vcfs = population_vcf,
			population_vcf_indices = population_vcf_index,
			gff = select_first([reference.gff]),
			runtime_attributes = default_runtime_attributes
	}

	call slivar_svpack_tsv {
		input:
			filtered_vcf = svpack_filter_annotated.svpack_vcf,
			pedigree = write_yaml_ped_phrank.pedigree,
			lof_lookup = slivar_data.lof_lookup,
			clinvar_lookup = slivar_data.clinvar_lookup,
			phrank_lookup = write_yaml_ped_phrank.phrank_lookup,
			runtime_attributes = default_runtime_attributes
	}

	output {
		IndexData filtered_small_variant_vcf = {"data": slivar_small_variant.filtered_vcf, "data_index": slivar_small_variant.filtered_vcf_index}
		IndexData compound_het_small_variant_vcf = {"data": slivar_small_variant.compound_het_vcf, "data_index": slivar_small_variant.compound_het_vcf_index}
		File filtered_small_variant_tsv = slivar_small_variant.filtered_tsv
		File compound_het_small_variant_tsv = slivar_small_variant.compound_het_tsv
		IndexData filtered_svpack_vcf = {"data": svpack_filter_annotated.svpack_vcf, "data_index": svpack_filter_annotated.svpack_vcf_index}
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

task write_yaml_ped_phrank {
	input {
		String cohort_id
		File cohort_json

		File hpo_terms
		File hpo_dag
		File hpo_annotations
		File ensembl_to_hgnc

		RuntimeAttributes runtime_attributes
	}

	Int disk_size = ceil((size(hpo_terms, "GB") + size(hpo_dag, "GB") + size(hpo_annotations, "GB") + size(ensembl_to_hgnc, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		parse_cohort.py \
			--cohort_json ~{cohort_json} \
			--write_cohort_yaml ~{cohort_id}.yml

		yaml2ped.py \
			~{cohort_id}.yml \
			~{cohort_id} \
			~{cohort_id}.ped

		calculate_phrank.py \
			~{hpo_terms} \
			~{hpo_dag} \
			~{hpo_annotations} \
			~{ensembl_to_hgnc} \
			~{cohort_id}.yml \
			~{cohort_id} \
			~{cohort_id}_phrank.tsv
	>>>

	output {
		File cohort_yaml = "~{cohort_id}.yml"
		File pedigree = "~{cohort_id}.ped"
		File phrank_lookup = "~{cohort_id}_phrank.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pyyaml@sha256:af6f0689a7412b1edf76bd4bf6434e7fa6a86192eebf19573e8618880d9c1dbb"
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

task slivar_small_variant {
	input {
		File vcf
		File vcf_index

		File pedigree

		File reference
		File reference_index

		File slivar_js
		File gnomad_af
		File hprc_af
		File gff

		File lof_lookup
		File clinvar_lookup
		File phrank_lookup

		RuntimeAttributes runtime_attributes
	}

	Float max_gnomad_af = 0.03
	Float max_hprc_af = 0.03
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
	Array[String] info_fields = [
		'gnomad_af',
		'hprc_af',
		'gnomad_nhomalt',
		'hprc_nhomalt',
		'gnomad_ac',
		'hprc_ac'
	]

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int threads = 8
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(gnomad_af, "GB") + size(hprc_af, "GB") + size(gff, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

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

		bcftools index \
			--threads ~{threads - 1} \
			~{vcf_basename}.norm.bcf

		# slivar has no version option
		slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+' 
		
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
			--vcf ~{vcf_basename}.norm.bcf \
			--ped ~{pedigree} \
		| bcftools csq \
			--local-csq \
			--samples - \
			--ncsq 40 \
			--gff-annot ~{gff} \
			--fasta-ref ~{reference} \
			- \
			--output-type z \
			--output ~{vcf_basename}.norm.slivar.vcf.gz

		bcftools index \
			--threads ~{threads - 1} \
			--tbi ~{vcf_basename}.norm.slivar.vcf.gz

		slivar \
			compound-hets \
			--skip ~{sep=',' skip_list} \
			--vcf ~{vcf_basename}.norm.slivar.vcf.gz \
			--sample-field comphet_side \
			--ped ~{pedigree} \
			--allow-non-trios \
		| add_comphet_phase.py \
		| bcftools view \
			--output-type z \
			--output ~{vcf_basename}.norm.slivar.compound_hets.vcf.gz

		bcftools index \
			--threads ~{threads - 1} \
			--tbi ~{vcf_basename}.norm.slivar.compound_hets.vcf.gz

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
			~{vcf_basename}.norm.slivar.vcf.gz \
		| sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
		> ~{vcf_basename}.norm.slivar.tsv

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
			~{vcf_basename}.norm.slivar.compound_hets.vcf.gz \
		| sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
		> ~{vcf_basename}.norm.slivar.compound_hets.tsv
	>>>

	output {
		File filtered_vcf = "~{vcf_basename}.norm.slivar.vcf.gz"
		File filtered_vcf_index = "~{vcf_basename}.norm.slivar.vcf.gz.tbi"
		File compound_het_vcf = "~{vcf_basename}.norm.slivar.compound_hets.vcf.gz"
		File compound_het_vcf_index = "~{vcf_basename}.norm.slivar.compound_hets.vcf.gz.tbi"
		File filtered_tsv = "~{vcf_basename}.norm.slivar.tsv"
		File compound_het_tsv = "~{vcf_basename}.norm.slivar.compound_hets.tsv"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/slivar@sha256:0a09289ccb760da310669906c675be02fd16b18bbedc971605a587275e34966c"
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

		bgzip --version

		bgzip ~{sv_vcf_basename}.svpack.vcf

		tabix --version

		tabix -p vcf ~{sv_vcf_basename}.svpack.vcf.gz
	>>>

	output {
		File svpack_vcf = "~{sv_vcf_basename}.svpack.vcf.gz"
		File svpack_vcf_index = "~{sv_vcf_basename}.svpack.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/svpack@sha256:a680421cb517e1fa4a3097838719a13a6bd655a5e6980ace1b03af9dd707dd75"
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
		docker: "~{runtime_attributes.container_registry}/slivar@sha256:0a09289ccb760da310669906c675be02fd16b18bbedc971605a587275e34966c"
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
