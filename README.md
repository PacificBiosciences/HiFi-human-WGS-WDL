<h1 align="center"><img width="300px" src="https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/raw/main/images/logo_wdl_workflows.svg"/></h1>

<h1 align="center">PacBio WGS Variant Pipeline</h1>

Workflow for analyzing human PacBio whole genome sequencing (WGS) data using [Workflow Description Language (WDL)](https://openwdl.org/).

- Docker images used by this workflow are defined in [the wdl-dockerfiles repo](https://github.com/PacificBiosciences/wdl-dockerfiles). Images are hosted in PacBio's [quay.io](https://quay.io/organization/pacbio).
- Common tasks that may be reused within or between workflows are defined in [the wdl-common repo](https://github.com/PacificBiosciences/wdl-common).

# Workflow

**Workflow entrypoint**: [workflows/main.wdl](workflows/main.wdl)

PacBio WGS Variant Pipeline performs read alignment, variant calling, and phasing. Joint-calling of small variants and structural variants for cohorts and optional variant filtering and annotation is also available for HiFi human WGS. The workflow can run using Azure, AWS, GCP, and HPC backends.

![PacBio WGS Variant Pipeline diagram](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/raw/main/images/main.graphviz.svg "PacBio WGS Variant Pipeline diagram")

## Setup

We recommend cloning the repo rather than downloading the release package.  Some tasks and workflows are pulled in from other repositories. Ensure you have initialized submodules following cloning by running `git submodule update --init --recursive`.

## Resource requirements

The workflow requires at minimum 64 cores and 256 GB of RAM. Ensure that the backend environment you're using has enough quota to run the workflow.

## Reference datasets and associated workflow files

Reference datasets are hosted publicly for use in the pipeline. For data locations, see the [backend-specific documentation](backends/) and template inputs files for each backend with paths to publicly hosted reference files filled out.

# Running the workflow

1. [Select a backend environment](#selecting-a-backend)
2. [Configure a workflow execution engine in the chosen environment](#configuring-a-workflow-engine)
3. [Fill out the inputs JSON file for your cohort](#filling-out-the-inputs-json)
4. [Run the workflow](#running-the-workflow-1)

## Selecting a backend

The workflow can be run on Azure, AWS, GCP, or HPC. Your choice of backend will largely be determined by the location of your data.

For backend-specific configuration, see the relevant documentation:

- [Azure](backends/azure)
- [AWS](backends/aws)
- [GCP](backends/gcp)
- [HPC](backends/hpc)

## Configuring a workflow engine and container runtime

An execution engine is required to run workflows. Two popular engines for running WDL-based workflows are [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/getting_started.html) and [`Cromwell`](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

Because workflow dependencies are containerized, a container runtime is required. This workflow has been tested with [Docker](https://docs.docker.com/get-docker/) and [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/) container runtimes.

See the [backend-specific documentation](backends) for details on setting up an engine.

| Engine | Azure | AWS | GCP | HPC |
| :- | :- | :- | :- | :- |
| [**miniwdl**](https://github.com/chanzuckerberg/miniwdl#scaling-up) | _Unsupported_ | Supported via the [Amazon Genomics CLI](https://aws.amazon.com/genomics-cli/) | _Unsupported_ | (SLURM only) Supported via the [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm) plugin |
| [**Cromwell**](https://cromwell.readthedocs.io/en/stable/backends/Backends/) | Supported via [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) | Supported via the [Amazon Genomics CLI](https://aws.amazon.com/genomics-cli/) | Supported via Google's [Pipelines API](https://cromwell.readthedocs.io/en/stable/backends/Google/) | Supported - [Configuration varies depending on HPC infrastructure](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/) |

## Filling out the inputs JSON

The input to a workflow run is defined in JSON format. Template input files with reference dataset information filled out are available for each backend:

- [Azure](backends/azure/inputs.azure.json)
- [AWS](backends/aws/inputs.aws.json)
- [GCP](backends/gcp/inputs.gcp.json)
- [HPC](backends/hpc/inputs.hpc.json)

Using the appropriate inputs template file, fill in the cohort and sample information (see [Workflow Inputs](#workflow-inputs) for more information on the input structure).

If using an HPC backend, you will need to download the reference bundle and replace the `<local_path_prefix>` in the input template file with the local path to the reference datasets on your HPC.

## Running the workflow

Run the workflow using the engine and backend that you have configured ([miniwdl](#run-directly-using-miniwdl), [Cromwell](#run-directly-using-cromwell)).

Note that the calls to `miniwdl` and `Cromwell` assume you are accessing the engine directly on the machine on which it has been deployed. Depending on the backend you have configured, you may be able to submit workflows using different methods (e.g. using trigger files in Azure, or using the Amazon Genomics CLI in AWS).

### Run directly using miniwdl

`miniwdl run workflows/main.wdl -i <input_file_path.json>`

### Run directly using Cromwell

`java -jar <cromwell_jar_path> run workflows/main.wdl -i <input_file_path.json>`

If Cromwell is running in server mode, the workflow can be submitted using cURL. Fill in the values of CROMWELL_URL and INPUTS_JSON below, then from the root of the repository, run:

```bash
# The base URL (and port, if applicable) of your Cromwell server
CROMWELL_URL=
# The path to your inputs JSON file
INPUTS_JSON=

(cd workflows && zip -r dependencies.zip humanwgs_structs.wdl  cohort_analysis/ sample_analysis/ tertiary_analysis/ wdl-common/)
curl -X "POST" \
  "${CROMWELL_URL}/api/workflows/v1" \
  -H "accept: application/json" \
  -H "Content-Type: multipart/form-data" \
  -F "workflowSource=@workflows/main.wdl" \
  -F "workflowInputs=@${INPUTS_JSON};type=application/json" \
  -F "workflowDependencies=@workflows/dependencies.zip;type=application/zip"
```

To specify [workflow options](https://cromwell.readthedocs.io/en/latest/wf_options/Overview/), add the following to the request (assuming your options file is a file called `options.json` located in the `pwd`): `-F "workflowOptions=@options.json;type=application/json"`.

# Workflow inputs

This section describes the inputs required for a run of the workflow. Typically, only the `humanwgs.cohort` and potentially [run/backend-specific sections](#other-inputs) will be filled out by the user for each run of the workflow. Input templates with reference file locations filled out are provided [for each backend](backends).

## [Cohort](workflows/humanwgs_structs.wdl)

A cohort can include one or more samples. Samples need not be related, but if you plan to run tertiary analysis, it is best to think of a cohort as a family of related samples. We have tested cohorts of up to 5 samples with 30x coverage.  Larger cohorts may encounter memory issues during joint calling.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | cohort_id | A unique name for the cohort; used to name outputs | |
| Array[[Sample](#sample)] | samples | The set of samples for the cohort. At least one sample must be defined. | |
| Array[String] | phenotypes | [Human Phenotype Ontology (HPO) phenotypes](https://hpo.jax.org/app/) associated with the cohort. If no particular phenotypes are desired, the root HPO term, `"HP:0000001"`, can be used. | |

### [Sample](workflows/humanwgs_structs.wdl)

Sample information for each sample in the workflow run.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | sample_id | A unique name for the sample; used to name outputs | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | movie_bams | The set of unaligned movie BAMs associated with this sample | |
| String? | sex | Sample sex | ["MALE", "FEMALE", `null`]. If the sex field is missing or `null`, sex will be set to unknown. Used to set the expected sex chromosome karyotype for TRGT and HiFiCNV. |
| Boolean | affected | Is this sample affected by the phenotype? | \[`true`, `false`\] |
| String? | father_id | Paternal `sample_id` | |
| String? | mother_id | Maternal `sample_id` | |

## [ReferenceData](workflows/humanwgs_structs.wdl)

Files associated with the reference genome.

These files are hosted publicly in each of the cloud backends; see `backends/${backend}/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | name | Reference name; used to name outputs (e.g., "GRCh38") | Note: The workflow currently only supports GRCh38 and provides GCA_000001405.15_GRCh38_no_alt_analysis_set. |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | fasta | Reference genome and index | |
| File | tandem_repeat_bed | Tandem repeat locations used by [pbsv](https://github.com/PacificBiosciences/pbsv) to normalize SV representation | |
| File | trgt_tandem_repeat_bed | Tandem repeat sites to be genotyped by [TRGT](https://github.com/PacificBiosciences/trgt) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | hificnv_exclude_bed | Compressed BED and index of regions to exclude from calling by [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV).  We recommend [cnv.excluded_regions.common_50.hg38.bed.gz](https://github.com/PacificBiosciences/HiFiCNV/blob/main/docs/aux_data.md). | |
| File | hificnv_expected_bed_male | BED of expected copy number for male karyotype for HiFiCNV | |
| File | hificnv_expected_bed_female | BED of expected copy number for female karyotype for HiFiCNV | |
| File? | gnomad_af | [gnomAD](https://gnomad.broadinstitute.org/) v3.1 allele frequences in [`slivar gnotate`](https://github.com/brentp/slivar/wiki/gnotate) format | required if `run_tertiary_analysis` is set to `true` |
| File? | hprc_af | Allele frequences in ~100 [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/) samples in `slivar gnotate` format | required if `run_tertiary_analysis` is set to `true` |
| File? | gff | [Ensembl](https://useast.ensembl.org/index.html) GFF3 reference annotation | required if `run_tertiary_analysis` is set to `true` |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)?] | population_vcfs | An array of structural variant population VCFs | required if `run_tertiary_analysis` is set to `true` |

## [SlivarData](workflows/humanwgs_structs.wdl)

Files associated with `slivar` annotation.  These are required if `run_tertiary_analysis` is set to `true`.

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
| Boolean? | run_tertiary_analysis | Run the optional tertiary analysis steps \[`false`\] | \[`true`, `false`\] |
| String | backend | Backend where the workflow will be executed | \["Azure", "AWS", "GCP", "HPC"\] |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'. | <ul><li>[Determining available zones in AWS](backends/aws/README.md#determining-available-zones)</li><li>[Determining available zones in GCP](backends/gcp/README.md#determining-available-zones)</li></ul> |
| String? | aws_spot_queue_arn | Queue ARN for the spot batch queue; required if backend is set to 'AWS' and `preemptible` is set to `true` | [Determining the AWS queue ARN](backends/aws/README.md#determining-the-aws-batch-queue-arn) |
| String? | aws_on_demand_queue_arn | Queue ARN for the on demand batch queue; required if backend is set to 'AWS' and `preemptible` is set to `false` | [Determining the AWS queue ARN](backends/aws/README.md#determining-the-aws-batch-queue-arn) |
| String? | container_registry | Container registry where workflow images are hosted. If left blank, [PacBio's public Quay.io registry](https://quay.io/organization/pacbio) will be used. | |
| Boolean | preemptible | If set to `true`, run tasks preemptibly where possible. On-demand VMs will be used only for tasks that run for >24 hours if the backend is set to GCP. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. | \[`true`, `false`\] |

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
| Array[File] | small_variant_roh_out | Output of [`bcftools roh`](https://samtools.github.io/bcftools/howtos/roh-calling.html) using `--AF-dflt 0.4` | |
| Array[File] | small_variant_roh_bed | Regions of homozygosity determiend by [`bcftools roh`](https://samtools.github.io/bcftools/howtos/roh-calling.html) using `--AF-dflt 0.4` | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_phased_small_variant_vcfs | Small variants called by DeepVariant and phased by [HiPhase](https://github.com/PacificBiosciences/HiPhase) (with index) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_phased_sv_vcfs | Structural variants called by [pbsv](https://github.com/PacificBiosciences/pbsv) and phased by HiPhase (with index) | |
| Array[File] | sample_hiphase_stats | Phase block summary statistics written by [HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#chromosome-summary-file---summary-file) | |
| Array[File] | sample_hiphase_blocks | Phase block list written by [HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#phase-block-file---blocks-file) | |
| Array[File] | sample_hiphase_haplotags | Per-read haplotag information, written by [HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#haplotag-file---haplotag-file) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | merged_haplotagged_bam | Aligned (by [pbmm2](https://github.com/PacificBiosciences/pbmm2)), haplotagged (by [HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#haplotagged-bam-files)) reads (with index) | |
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
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_small_variant_vcf | Small variants called by [DeepVariant](https://github.com/google/deepvariant), joint-called by [GLnexus](https://github.com/dnanexus-rnd/GLnexus), and phased by [HiPhase](https://github.com/PacificBiosciences/HiPhase) (with index) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_sv_vcf | Structural variants joint-called by [pbsv](https://github.com/PacificBiosciences/pbsv) and phased by HiPhase (with index) | |
| File? | cohort_hiphase_stats | Phase block summary statistics written by [HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#chromosome-summary-file---summary-file) | |
| File? | cohort_hiphase_blocks | Phase block list written by [HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md#phase-block-file---blocks-file) | |

## Tertiary analysis

These files will be output for each run of the workflow if `run_tertiary_analysis` is set to `true`. The files that are being annotated will depend on whether the number of samples is equal to or greater than one:
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

Docker images definitions used by this workflow can be found in [the wdl-dockerfiles repository](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a). Images are hosted in PacBio's [quay.io](https://quay.io/organization/pacbio). Docker images used in the workflow are pegged to specific versions by referring to their digests rather than tags.

The Docker image used by a particular step of the workflow can be identified by looking at the `docker` key in the `runtime` block for the given task. Images can be referenced in the following table by looking for the name after the final `/` character and before the `@sha256:...`. For example, the image referred to here is "align_hifiasm":
> ~{runtime_attributes.container_registry}/**align_hifiasm**@sha256:3968cb<...>b01f80fe

| Image | Major tool versions | Links |
| :- | :- | :- |
| bcftools | <ul><li>[bcftools 1.14](https://github.com/samtools/bcftools/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/bcftools) |
| deepvariant | User-defined; default is version [1.5.0](https://github.com/google/deepvariant/releases/tag/v1.5.0) | [DeepVariant GitHub](https://github.com/google/deepvariant) |
| glnexus | <ul><li>[glnexus v1.4.3](https://github.com/dnanexus-rnd/GLnexus/releases/tag/v1.4.3)</li></ul> | [GLnexus GitHub](https://github.com/dnanexus-rnd/GLnexus) |
| hificnv | <ul><li>[HiFiCNV v0.1.7](https://github.com/PacificBiosciences/HiFiCNV/releases/tag/v0.1.7)</li><li>[bcftools 1.16](https://github.com/samtools/bcftools/releases/tag/1.16)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/0b0fbe939648087e9fdea4497ae08dc76538ebf0/docker/hificnv) |
| hiphase | <ul><li>[HiPhase 1.0.0](https://github.com/PacificBiosciences/HiPhase/releases/tag/v1.0.0)</li><li>[samtools 1.18](https://github.com/samtools/samtools/releases/tag/1.18)</li><li>[bcftools 1.18](https://github.com/samtools/bcftools/releases/tag/1.18)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/d26db6204409dfeff56e169cdba0cc14bc272f15/docker/hiphase) |
| htslib | <ul><li>[htslib 1.14](https://github.com/samtools/htslib/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/htslib) |
| mosdepth | <ul><li>[mosdepth 0.2.9](https://github.com/brentp/mosdepth/releases/tag/v0.2.9)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/mosdepth) |
| paraphase | <ul><li>[minimap2 2.17](https://github.com/lh3/minimap2/releases/tag/v2.17)</li><li>[samtools 1.14](https://github.com/samtools/samtools/releases/tag/1.14)</li><li>[paraphase 2.2.3](https://github.com/PacificBiosciences/paraphase/releases/tag/v2.2.3)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/paraphase) |
| pb-cpg-tools | <ul><li>[pb-CpG-tools v2.3.2](https://github.com/PacificBiosciences/pb-CpG-tools/releases/tag/v2.3.2)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/7481837d3b0f539adf4f64209a65cf28eebf3dba/docker/pb-cpg-tools) |
| pbmm2 | <ul><li>[pbmm2 1.10.0](https://github.com/PacificBiosciences/pbmm2/releases/tag/v1.10.0)</li><li>[datamash 1.1.0](https://ftp.gnu.org/gnu/datamash/)</li><li>[pysam 0.16.0.1](https://github.com/pysam-developers/pysam/releases/tag/v0.16.0.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/pbmm2) |
| pbsv | <ul><li>[pbsv 2.9.0](https://github.com/PacificBiosciences/pbsv/releases/tag/v2.9.0)</li><li>[htslib 1.14](https://github.com/samtools/htslib/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/f9e33a757e6d8cb15696ac930a2efd0fd7a885d8/docker/pbsv) |
| pyyaml | <ul><li>[pyyaml 5.3.1](https://github.com/yaml/pyyaml/releases/tag/5.3.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/f72e862bca2f209b9909e6043ef0197975762f27/docker/pyyaml) |
| samtools | <ul><li>[samtools 1.14](https://github.com/samtools/samtools/releases/tag/1.14)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/samtools) |
| slivar | <ul><li>[slivar 0.2.2](https://github.com/brentp/slivar/releases/tag/v0.2.2)</li><li>[bcftools 1.14](https://github.com/samtools/bcftools/releases/tag/1.14)</li><li>[vcfpy 0.13.3](https://github.com/bihealth/vcfpy/releases/tag/v0.13.3)</li><li>[pysam 0.19.1](https://github.com/pysam-developers/pysam/releases/tag/v0.19.1)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3560fcc5a84e044067cea9c9a7669cfc2659178e/docker/slivar) |
| svpack | <ul><li>[svpack 36180ae6](https://github.com/PacificBiosciences/svpack/tree/a82598ebc4013bf32e70295b83b380ada6302c4a)</li><li>[htslib 1.18](https://github.com/samtools/htslib/releases/tag/1.18)</li><li>[pysam 0.21.0](https://github.com/pysam-developers/pysam/releases/tag/v0.21.0)</li> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/8edbc516abc0ff43ac279b48018003923721b054/docker/svpack) |
| trgt | <ul><li>[trgt 0.5.0](https://github.com/PacificBiosciences/trgt/releases/tag/v0.5.0)</li><li>[samtools 1.18](https://github.com/samtools/samtools/releases/tag/1.18)</li><li>[bcftools 1.18](https://github.com/samtools/bcftools/releases/tag/1.18)</li><li>[pysam 0.21.0](https://github.com/pysam-developers/pysam/releases/tag/v0.21.0)</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/d2a45e0213ac3fa631a51a48757c442d3ed550b6/docker/trgt) |

---

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
