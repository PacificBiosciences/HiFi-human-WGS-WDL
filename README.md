# wdl-humanwgs

Workflow for analyzing human PacBio whole genome sequencing (WGS) data using [Workflow Description Language (WDL)](https://openwdl.org/).

- For the snakemake version of these workflows, see [here](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake).

- Docker images used by these workflows are defined in [the wdl-dockerfiles repo](https://github.com/PacificBiosciences/wdl-dockerfiles). Images are hosted in PacBio's [quay.io](https://quay.io/organization/pacbio).

- Common tasks that may be reused within or between workflows are defined in [the wdl-common repo](https://github.com/PacificBiosciences/wdl-common). This repository is a submodule of the wdl-humanwgs repo; ensure the submodule has been initialized prior to attempting to run the workflow  (`git submodule update --init --recursive`).

# Workflow

This human WGS workflow performs alignment, variant calling, and variant phasing, as well as joint-calling small variants and structural variants for families, and optional variant filtering and annotation. The workflow can run using Azure, AWS, GCP, and HPC backends.

![Human WGS workflow diagram](workflows/main.graphviz.svg "Human WGS workflow diagram")

# Running the workflow

**Workflow entrypoint**: [workflows/main.wdl](workflows/main.wdl)

Some tasks and workflows are pulled in from other repositories. Make sure you have initialized submodules following cloning by running `git submodule update --init --recursive`.

## Backend environments

The workflow can be run on Azure, AWS, GCP, or HPC. For backend-specific configuration, see the relevant documentation:

- [Azure](backends/azure)
- [AWS](backends/aws)
- [GCP](backends/gcp)
- [HPC](backends/hpc)

## Resource requirements

The workflow requires at minimum 64 cores and 48 GB of RAM. Ensure that the backend environment you're using has enough quota to run the workflow.

## Workflow engines

Two popular engines for running WDL-based workflows are [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/getting_started.html) and [`Cromwell`](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

The workflow engine that you choose will depend on where your data is located.

| Engine | Azure | AWS | GCP | HPC |
| :- | :- | :- | :- | :- |
| [**miniwdl**](https://github.com/chanzuckerberg/miniwdl#scaling-up) | _Unsupported_ | Supported via the [Amazon Genomics CLI](https://aws.amazon.com/genomics-cli/) | _Unsupported_ | (SLURM only) Supported via the [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm) plugin |
| [**Cromwell**](https://cromwell.readthedocs.io/en/stable/backends/Backends/) | Supported via [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) | Supported via the [Amazon Genomics CLI](https://aws.amazon.com/genomics-cli/) | Supported via Google's [Pipelines API](https://cromwell.readthedocs.io/en/stable/backends/Google/) | Supported - [Configuration varies depending on HPC infrastructure](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/) |

## Setting up and running the workflow

1. Install and configure the workflow execution engine of your choice following the documentation for the backend environment where your data is located. See the [backend environments](#backend-environments) section for more information on configuring engines and quotas.

2. Select the input template file ([Azure](backends/azure/inputs.azure.json), [AWS](backends/aws/inputs.aws.json), [GCP](backends/gcp/inputs.gcp.json), [HPC](backends/hpc/inputs.hpc.json)) that matches the backend environment where you will be executing the workflows. These files have reference dataset information prefilled. If using an HPC backend, you will need to download the reference bundle and replace the `<local_path_prefix>` in the input template file with the local path to the reference datasets on your HPC.

3. Fill in the cohort and sample information (see [Workflow Inputs](#workflow-inputs) for more information on the input structure).

4. Run the workflow using the engine and backend of choice ([miniwdl](#running-using-miniwdl), [Cromwell](#running-using-cromwell), [Workbench](#running-using-workbench)).

### Run using miniwdl

`miniwdl run workflows/main.wdl -i <input_file_path.json>`

### Run using Cromwell

`java -jar <cromwell_jar_path> run workflows/main.wdl -i <input_file_path.json>`

### Run using Workbench

```bash
dnastack alpha workbench runs submit \
	--workflow-params '@<input_file_path.json>' \
	--url <internalId>
```

See the next section for details on configuring engines and submitting workflow runs using Workbench.

## Running and monitoring workflows using Workbench

Rather than running a workflow directly using an engine, engines can be configured using [Workbench](https://workbench.dnastack.com/), a software suite that enables users to easily configure, run, and monitor workflow runs across backend environments. Engines are set up in the desired backend and then registered with Workbench, following which runs may be submitted and monitored either directly in-browser or using the command-line interface (CLI).

Note that a license is required to submit runs using Workbench.

### Configuring an engine

Workflows submitted via Workbench must be submitted to an execution engine. See [the Workbench documentation](https://docs.dnastack.com/docs/workbench-settings) for details about setting up an engine in the backend of your choice. Backend-specific resources and default configurations that may be required as part of engine setup may also be found in the [backends](backends) directory.

Once the engine has been set up, follow the documentation for registering the engine with Workbench, following which you will be able to submit workflow runs using Workbench either [via the browser](https://docs.dnastack.com/docs/accessing-the-workbench-gui) or [via the CLI](#running-the-workflow-using-the-workbench-cli).

### Installing and configuring the DNAstack CLI

1. Install the DNAstack CLI

`python3 -m pip install --user dnastack-client-library`

Confirm that the CLI is installed and available by running `dnastack --version`.

2. Authenticate using the CLI

`dnastack auth login`

3. Configure the CLI to use workbench

`dnastack use workbench.dnastack.com`

You can now use the DNAstack CLI to interact with Workbench. See [the CLI documentation](https://docs.dnastack.com/docs/runs-submit) for details on available commands.

### Running the workflow using the Workbench CLI

Note that these steps assume you have already [set up and registered an engine in Workbench](https://docs.dnastack.com/docs/workbench-settings).

1. Register the workflow on Workbench

From the root of this repository, run:

```bash
dnastack alpha workbench workflows create \
	--name "PacBio HumanWGS" \
	--description =@README.md \
	workflows/main.wdl
```

Note the `internalId` field of the returned JSON. This will be used as the `--url` value when submitting workflow runs.

2. Fill out the inputs.json file

The inputs you use will depend on the backend where you've set up your engine. Template files for [AWS](backends/aws/inputs.aws.json), [Azure](backends/azure/inputs.azure.json), and [GCP](backends/gcp/inputs.gcp.json) with reference file locations already filled out can be used to get started; see [the workflow inputs section](#workflow-inputs) for information on the other fields you will need to fill out.

3. Submit the workflow to Workbench

In the following command, replace `<input_file_path.json>` with the path to your filled out inputs file, and `<internalId>` with the ID you noted in step 1. If no engine is provided, the default engine you have configured will be used.

```bash
dnastack alpha workbench runs submit \
	--workflow-params '@<input_file_path.json>' \
	--url <internalId> \
	[--tags <key=value>] \
	[--engine <engineId>]
```

# Reference datasets and associated workflow files

Reference datasets are hosted publicly for use in the pipeline. For data locations, the [backend-specific documentation](backends/) and template inputs files for each backend with paths to publicly hosted reference files filled out.

# Workflow inputs

## [Cohort](workflows/humanwgs_structs.wdl)

A cohort can include one or more samples. Samples need not be related, but if you plan to run tertiary analysis, it is best to think of a cohort as a family of related samples. We have tested cohorts of up to 5 samples with 30x coverage.  Larger cohorts may encounter memory issues during joint calling.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | cohort_id | A unique name for the cohort; used to name outputs | |
| Array[[Sample](#sample)] | samples | The set of samples for the cohort. At least one sample must be defined. | |
| Array[String] | phenotypes | [Human Phenotype Ontology (HPO) phenotypes](https://hpo.jax.org/app/) associated with the cohort | |

### [Sample](workflows/humanwgs_structs.wdl)

Sample information for each sample in the workflow run.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | sample_id | A unique name for the sample; used to name outputs | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | movie_bams | The set of unaligned movie BAMs associated with this sample | |
| String? | sex | Sample sex | ["MALE", "FEMALE", `null`]. If the sex field is missing or `null`, sex will be set to unknown. Used to set the expected sex chromosome karyotype for TRGT and HiFiCNV. |
| Boolean | affected | Is this sample affected by the phenotype? | \[true, false\] |
| String? | father_id | Paternal `sample_id` | |
| String? | mother_id | Maternal `sample_id` | |

## [ReferenceData](workflows/humanwgs_structs.wdl)

Files associated with the reference genome.

These files are hosted publicly in each of the cloud backends; see `backends/${backend}/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | name | Reference name; used to name outputs (e.g., "GRCh38") | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | fasta | Reference genome and index | |
| Array[String] | chromosomes | Chromosomes to phase, typically `chr{1..22} chr{X,Y}` | |
| File | chromosome_lengths | Reference chromosome lengths | |
| File | tandem_repeat_bed | Tandem repeat locations used by [pbsv](https://github.com/PacificBiosciences/pbsv) to normalize SV representation | |
| File | trgt_tandem_repeat_bed | Tandem repeat sites to be genotyped by [TRGT](https://github.com/PacificBiosciences/trgt) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | hificnv_exclude_bed | Compressed BED and index of regions to exclude from calling by [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV).  We recommend [cnv.excluded_regions.common_50.hg38.bed.gz](https://github.com/PacificBiosciences/HiFiCNV/blob/main/docs/aux_data.md). | |
| File | hificnv_expected_bed_male | BED of expected copy number for male karyotype for HiFiCNV | |
| File | hificnv_expected_bed_female | BED of expected copy number for female karyotype for HiFiCNV | |
| File | gnomad_af | [gnomAD](https://gnomad.broadinstitute.org/) v3.1 allele frequences in [`slivar gnotate`](https://github.com/brentp/slivar/wiki/gnotate) format | |
| File | hprc_af | Allele frequences in ~100 [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/) samples in `slivar gnotate` format | |
| File | gff | [Ensembl](https://useast.ensembl.org/index.html) GFF3 reference annotation | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | population_vcfs | An array of structural variant population VCFs | |

## [SlivarData](workflows/humanwgs_structs.wdl)

Files associated with `slivar` annotation.

These files are hosted publicly in each of the cloud backends; see `backends/${backend}/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| File | slivar_js | Additional javascript functions for slivar | |
| File | hpo_terms | [HPO](https://hpo.jax.org/app/) annotation lookups | |
| File | hpo_dag | HPO annotation lookups | |
| File | hpo_annotations | HPO annotation lookups | |
| File | ensembl_to_hgnc | Ensembl to HGNC gene mapping | |
| File | lof_lookup | Loss-of-function scores per gene | |
| File | clinvar_lookup | ClinVar annotations per gene | |

## Other inputs

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String? | deepvariant_version | Version of deepvariant to use \["1.5.0"\] | |
| [DeepVariantModel](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | deepvariant_model | Optional alternate DeepVariant model file to use | |
| Int? | pbsv_call_mem_gb | Optionally set RAM (GB) for pbsv_call during cohort analysis | |
| Int? | glnexus_mem_gb | Optionally set RAM (GB) for GLnexus during cohort analysis | |
| Boolean? | run_tertiary_analysis | Run the optional tertiary analysis steps \[true\] | |
| String | backend | Backend where the workflow will be executed | \["Azure", "AWS", "GCP", "HPC"\] |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'. | <ul><li>[Determining available zones in AWS](backends/aws/README.md#determining-available-zones)</li><li>[Determining available zones in GCP](backends/gcp/README.md#determining-available-zones)</li></ul> |
| String? | aws_spot_queue_arn | Queue ARN for the spot batch queue; required if backend is set to 'AWS' and `preemptible` is set to `true` | [Determining the AWS queue ARN](backends/aws/README.md#determining-the-aws-batch-queue-arn) |
| String? | aws_on_demand_queue_arn | Queue ARN for the on demand batch queue; required if backend is set to 'AWS' and `preemptible` is set to `false` | [Determining the AWS queue ARN](backends/aws/README.md#determining-the-aws-batch-queue-arn) |
| Boolean | preemptible | If set to `true`, run tasks preemptibly where possible. On-demand VMs will be used only for tasks that run for >24 hours if the backend is set to GCP. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. | \[true, false\] |

# Workflow outputs

## Sample analysis

These files will be output for each sample defined in the cohort.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Array[Array[File]] | bam_stats | TSV of length and quality for each read (per input BAM) | |
| Array[Array[File]] | read_length_summary | For each input BAM, read length distribution (per input BAM) | |
| Array[Array[File]] | read_quality_summary | For each input BAM, read quality distribution (per input BAM) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | small_variant_gvcfs | Small variants (SNPs and INDELs < 50bp) gVCFs called by [DeepVariant](https://github.com/google/deepvariant) (with index) | |
| Array[File] | small_variant_vcf_stats | [`bcftools stats`](https://samtools.github.io/bcftools/bcftools.html#stats) summary statistics for small variants | |
| Array[File] | small_variant_roh_bed | Regions of homozygosity determiend by [`bcftools roh`](https://samtools.github.io/bcftools/howtos/roh-calling.html) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_sv_vcfs | Structural variants called by [pbsv](https://github.com/PacificBiosciences/pbsv) (with index) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_phased_small_variant_vcfs | Small variants called by DeepVariant and phased by [WhatsHap](https://whatshap.readthedocs.io/en/latest/) (with index) | |
| Array[File] | sample_whatshap_stats_tsvs | Phase block statistics written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats) | |
| Array[File] | sample_whatshap_stats_gtfs | Phase block GTF written by `whatshap stats` | |
| Array[File] | sample_whatshap_stats_blocklists | Haplotype block list written by `whatshap stats` | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | merged_haplotagged_bam | Aligned (by [pbmm2](https://github.com/PacificBiosciences/pbmm2)), haplotagged (by [`whatshap haplotag`](https://whatshap.readthedocs.io/en/latest/guide.html#visualizing-phasing-results)) reads (with index) | |
| Array[File] | haplotagged_bam_mosdepth_summary | [mosdepth](https://github.com/brentp/mosdepth) summary of median depths per chromosome | |
| Array[File] | haplotagged_bam_mosdepth_region_bed | mosdepthhttps://github.com/brentp/mosdepth BED of median coverage depth per 500 bp window | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | trgt_repeat_vcf | Tandem repeat genotypes from [TRGT](https://github.com/PacificBiosciences/trgt/blob/main/docs/vcf_files.md) (with index) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | trgt_spanning_reads | Fragments of HiFi reads spanning loci genotyped by TRGT (with index) | |
| Array[File] | trgt_dropouts | Regions with insufficient coverage to genotype by TRGT | |
| Array[Array[File]] | cpg_pileup_beds | 5mCpG site methylation probability pileups generated by [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools#output-files) | |
| Array[Array[File]] | cpg_pileup_bigwigs | 5mCpG site methylation probability pileups generated by pb-CpG-tools | |
| Array[File] | paraphase_output | Output generated by [Paraphase](https://github.com/PacificBiosciences/paraphase) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | paraphase_realigned_bam | Realigned BAM for selected medically relevant genes in segmental duplications (with index), generated by Paraphase | |
| Array[Array[File]] | paraphase_vcfs | Phased Variant calls for selected medically relevant genes in segmental duplications, generated by Paraphase | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | hificnv_vcfs | VCF output containing copy number variant calls for the sample from [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) | |
| Array[File] | hificnv_copynum_bedgraphs | Copy number values calculated for each region | |
| Array[File] | hificnv_depth_bws | Bigwig file containing the depth measurements from HiFiCNV | |
| Array[File] | hificnv_maf_bws | Bigwig file containing the minor allele frequency measurements from DeepVariant, generated by HiFiCNV | |

## Cohort analysis

These files will be output if the cohort includes more than one sample.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_sv_vcf | Structural variants joint-called by [pbsv](https://github.com/PacificBiosciences/pbsv) (with index) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_phased_joint_called_vcf | Small variants called by [DeepVariant](https://github.com/google/deepvariant), joint-called by [GLnexus](https://github.com/dnanexus-rnd/GLnexus), and phased by [WhatsHap](https://whatshap.readthedocs.io/en/latest/) (with index) | |
| File? | cohort_whatshap_stats_tsvs | Phase block statistics written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats)  | |
| File? | cohort_whatshap_stats_gtfs | Phase block GTF written by `whatshap stats` | |
| File? | cohort_whatshap_stats_blocklists | Haplotype block list written by `whatshap stats` | |

## Tertiary analysis

These files will be output for each run of the workflow if `run_tertiary_analysis` is set to `true` (this is the default). The files that are being annotated will depend on whether the number of samples is equal to or greater than one:
- If the number of samples is equal to one, the files being annotated in this step are the sample small variant VCF and SV VCF.
- If the number of samples is greater than one, the files being annotated in this step are the phased, joint-called small variant VCF and the cohort SV VCF.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | filtered_small_variant_vcf | Small variant calls that are filtered based on population frequency and annotated with cohort information, population frequency, gene, functional impact, etc., using [slivar](https://github.com/brentp/slivar) and [`bcftools csq`](https://samtools.github.io/bcftools/howtos/csq-calling.html) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | compound_het_small_variant_vcf | Compound heterozygotes annotated with cohort information, population frequency, gene, functional impact, etc., using slivar and `bcftools csq` | |
| File? | filtered_small_variant_tsv | Filtered VCFs are reformatted as a human-readable TSV by [`slivar tsv`](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF) | |
| File? | compound_het_small_variant_tsv | Filtered VCFs are reformatted as a human-readable TSV by `slivar tsv` | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | filtered_svpack_vcf | Structural variant calls that are filtered based on population frequency and annotated with cohort information, population frequency, gene, functional impact, etc., using [svpack](https://github.com/PacificBiosciences/svpack) | |
| File? | filtered_svpack_tsv | Filtered VCFs are reformatted as a human-readable TSV by `slivar tsv` | |

# Tool versions and Docker images

Docker images definitions used by the human WGS workflow can be found in [the wdl-dockerfiles repository](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a). Images are hosted in PacBio's [quay.io](https://quay.io/organization/pacbio). Docker images used in the workflow are pegged to specific versions by referring to their digests rather than tags.

The Docker image used by a particular step of the workflow can be identified by looking at the `docker` key in the `runtime` block for the given task. Images can be referenced in the following table by looking for the name after the final `/` character and before the `@sha256:...`. For example, the image referred to here is "align_hifiasm":
> ~{runtime_attributes.container_registry}/**align_hifiasm**@sha256:3968cb<...>b01f80fe

| Image | Major tool versions | Links |
| :- | :- | :- |
| bcftools | <ul><li>[bcftools 1.14](https://github.com/samtools/bcftools/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/bcftools) |
| deepvariant | User-defined; default is version [1.5.0](https://github.com/google/deepvariant/releases/tag/v1.5.0) | [DeepVariant GitHub](https://github.com/google/deepvariant) |
| glnexus | <ul><li>[glnexus v1.4.1](https://github.com/dnanexus-rnd/GLnexus/releases/tag/v1.4.1)</li></ul> | [GLnexus GitHub](https://github.com/dnanexus-rnd/GLnexus) |
| hificnv | <ul><li>[HiFiCNV v0.1.6](https://github.com/PacificBiosciences/HiFiCNV/releases/tag/v0.1.6)</li><li>[bcftools 1.16](https://github.com/samtools/bcftools/releases/tag/1.16)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/main/docker/hificnv) |
| htslib | <ul><li>[htslib 1.14](https://github.com/samtools/htslib/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/htslib) |
| mosdepth | <ul><li>[mosdepth 0.2.9](https://github.com/brentp/mosdepth/releases/tag/v0.2.9)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/mosdepth) |
| paraphase | <ul><li>[minimap2 2.17](https://github.com/lh3/minimap2/releases/tag/v2.17)</li><li>[samtools 1.14](https://github.com/samtools/samtools/releases/tag/1.14)</li><li>[paraphase 2.2.2](https://github.com/PacificBiosciences/paraphase/releases/tag/v2.2.2)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/ba57dca9958e936f1894065ab8e07263387bb476/docker/paraphase) |
| parse-cohort | <ul><li>python 3.8.10; custom scripts</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/parse-cohort) |
| pb-cpg-tools | <ul><li>[pb-CpG-tools v2.3.1](https://github.com/PacificBiosciences/pb-CpG-tools/releases/tag/v2.3.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/bb65fcbe05c3a42e6b37364c883471db47ef8de3/docker/pb-cpg-tools) |
| pbmm2 | <ul><li>[pbmm2 1.10.0](https://github.com/PacificBiosciences/pbmm2/releases/tag/v1.10.0)</li><li>[datamash 1.1.0](https://ftp.gnu.org/gnu/datamash/)</li><li>[pysam 0.16.0.1](https://github.com/pysam-developers/pysam/releases/tag/v0.16.0.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/pbmm2) |
| pbsv | <ul><li>[pbsv 2.9.0](https://github.com/PacificBiosciences/pbsv/releases/tag/v2.9.0)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/pbsv) |
| pyyaml | <ul><li>[pyyaml 5.3.1](https://github.com/yaml/pyyaml/releases/tag/5.3.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/pyyaml) |
| samtools | <ul><li>[samtools 1.14](https://github.com/samtools/samtools/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/samtools) |
| slivar | <ul><li>[slivar 0.2.2](https://github.com/brentp/slivar/releases/tag/v0.2.2)</li><li>[bcftools 1.14](https://github.com/samtools/bcftools/releases/tag/1.14)</li><li>[vcfpy 0.13.3](https://github.com/bihealth/vcfpy/releases/tag/v0.13.3)</li><li>[pysam 0.19.1](https://github.com/pysam-developers/pysam/releases/tag/v0.19.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/slivar) |
| svpack | <ul><li>[svpack 36180ae6](https://github.com/PacificBiosciences/svpack/tree/a82598ebc4013bf32e70295b83b380ada6302c4a)</li><li>[pysam 0.16.0.1](https://github.com/pysam-developers/pysam/releases/tag/v0.16.0.1)</li> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/8598b47771995e44799575b72bde477c21bfc213/docker/svpack) |
| trgt | <ul><li>[trgt 0.4.0](https://github.com/PacificBiosciences/trgt/releases/tag/v0.4.0)</li><li>[samtools 1.16.1](https://github.com/samtools/samtools/releases/tag/1.16.1)</li><li>[bcftools 1.16](https://github.com/samtools/bcftools/releases/tag/1.16)</li><li>[pysam 0.16.0.1](https://github.com/pysam-developers/pysam/releases/tag/v0.16.0.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/trgt) |
| whatshap | <ul><li>[whatshap 1.4](https://github.com/whatshap/whatshap/releases/tag/v1.4)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a/docker/whatshap) |
