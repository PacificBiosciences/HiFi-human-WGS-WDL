# wdl-humanwgs

Workflow for analyzing human PacBio whole genome sequencing (WGS) data using [Workflow Description Language (WDL)](https://openwdl.org/).

- For the snakemake version of these workflows, see [here](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake).

- Docker images used by these workflows are defined [here](https://github.com/PacificBiosciences/wdl-dockerfiles).

- Common tasks that may be reused within or between workflows are defined [here](https://github.com/PacificBiosciences/wdl-common).


# Workflow

The human WGS workflow performs read alignment, small and structural variant calling, variant phasing, and optional _de novo_ assembly.

**Workflow entrypoint**: [workflows/main.wdl](workflows/main.wdl)

- [Blank input template file](workflows/inputs.json)
- [Azure-based inputs](workflows/inputs.azure.json)
- [AWS-based inputs](workflows/inputs.aws.json)
- [GCP-based inputs]((workflows/inputs.gcp.json))

![Human WGS workflow diagram](workflows/main.graphviz.svg "Human WGS workflow diagram")


# Reference datasets and associated workflow files

Reference datasets are hosted publicly for use in the pipeline. For data locations, see `workflows/inputs.${backend}.json`.


## Reference data hosted in Azure

To use Azure reference data, add the following line to your `containers-to-mount` file in your Cromwell on Azure installation ([more info here](https://github.com/microsoft/CromwellOnAzure/blob/develop/docs/troubleshooting-guide.md#use-input-data-files-from-an-existing-azure-storage-account-that-my-lab-or-team-is-currently-using)):

`https://datasetpbrarediseases.blob.core.windows.net/dataset?si=public&spr=https&sv=2021-06-08&sr=c&sig=o6OkcqWWlGcGOOr8I8gCA%2BJwlpA%2FYsRz0DMB8CCtCJk%3D`

The [Azure input file template](workflows/inputs.azure.json) has paths to the reference files in this blob storage prefilled.


## Reference data hosted in AWS

AWS reference data is hosted in the `us-west-2` region  in the bucket `s3://dnastack-resources`.

To use AWS reference data, add the following line to the data section of your [`agc-project.yaml`](https://aws.github.io/amazon-genomics-cli/docs/concepts/projects/):

```
data:
  - location: s3://dnastack-resources
    readOnly: tue
```

The [AWS input file template](workflows/inputs.aws.json) has paths to the reference files in the blob storage prefilled.


## Reference data hosted in GCP

<TODO>


# Workflow inputs

## [Cohort](workflows/humanwgs_structs.wdl)

A cohort can include one or more samples. Samples need not be related.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | cohort_id | A unique name for the cohort; used to name outputs | |
| Array[[Sample](#sample)] | samples | The set of samples for the cohort. At least one sample must be defined. | |
| Array[String] | phenotypes | [HPO phenotypes](https://hpo.jax.org/app/) associated with the cohort | |
| Boolean | run_de_novo_assembly_trio | Run trio-based _de novo_ assembly. | Cohort must contain a valid trio (child and both parents present in the cohort) |


### [Sample](workflows/humanwgs_structs.wdl)

Sample information for each sample in the workflow run.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | sample_id | A unique name for the sample; used to name outputs | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | movie_bams | The set of movie bams associated with this sample | |
| String | sex | Sample sex | ["MALE", "FEMALE"] |
| Boolean | affected | The affected status for the sample | [true, false] |
| String? | father_id | Paternal `sample_id` | |
| String? | mother_id | Maternal `sample_id` | |
| Boolean | run_de_novo_assembly | If true, run single-sample _de novo_ assembly for this sample | [true, false] |


## [ReferenceData](workflows/humanwgs_structs.wdl)

Files associated with the reference genome.

These files are hosted publicly in each of the cloud backends; see `workflows/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | name | Reference name; used to name outputs | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | fasta | Reference genome and index to align reads to | |
| Array[String] | chromosomes | Chromosomes to phase during WhatsHap phasing | |
| File | chromosome_lengths | File specifying the lengths of each of the reference chromosomes | |
| File | tandem_repeat_bed | Tandem repeat locations in the reference genome | |
| File | trgt_tandem_repeat_bed | Repeat bed used by [`trgt`](https://github.com/PacificBiosciences/trgt) to output spanning reads and a repeat VCF | |
| File | gnomad_af | gnomAD allele frequencies; used for annotate the small variant VCF | |
| File | hprc_af | Allele frequences from the Human Pangenome Reference Consortium (HPRC); used to annotate the small variant VCF | |
| File | gff | GFF3 annotation file | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | population_vcfs | Population calls in VCF format; used to annotate the VCFs | |


## [SlivarData](workflows/humanwgs_structs.wdl)

Files associated with `slivar` annotation.

These files are hosted publicly in each of the cloud backends; see `workflows/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| File | slivar_js | Additional javascript functions for slivar | |
| File | hpo_terms | HPO annotation lookups | |
| File | hpo_dag | HPO annotation lookups | |
| File | hpo_annotations | HPO annotation lookups | |
| File | ensembl_to_hgnc | Ensembl to HGNC gene mapping | |
| File | lof_lookup | LoF lookup files for slivar annotations | |
| File | clinvar_lookup | ClinVar lookup files for slivar annotation | |


## Other inputs

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String? | deepvariant_version | Version of deepvariant to use [1.4.0] | |
| [DeepVariantModel](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | deepvariant_model | Optonal alternate DeepVariant model file to use | |
| String | backend | Backend where the workflow will be executed | ["Azure", "AWS", "GCP", "slurm"] |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'. | [Determining available zones in AWS and GCP](#determining-available-zones-in-aws-and-gcp). |
| String? | aws_spot_queue_arn | Queue ARN for the spot batch queue; required if backend is set to 'AWS' and `preemptible` is set to `true` | [Determining the AWS queue ARN](#determining-the-aws-batch-queue-arn) |
| String? | aws_on_demand_queue_arn | Queue ARN for the on demand batch queue; required if backend is set to 'AWS' and `preemptible` is set to `false` | [Determining the AWS queue ARN](#determining-the-aws-batch-queue-arn) |
| String? | slurm_partition_default | Default slurm partition; required if backend is set to 'slurm' |
| String? | slurm_partition_gpu | GPU slurm partition; optional if backend is set to 'slurm' |
| Boolean | preemptible | If set to `true`, run tasks preemptibly where possible. On-demand VMs will be used only for tasks that run for >24 hours if the backend is set to GCP. If set to `false`, on-demand VMs will be used for every task. | [true, false] |


### Determining available zones in AWS and GCP

#### AWS

To determine available zones in AWS, look for the ZoneName attributes output by the following command:

```bash
aws ec2 describe-availability-zones --region <region>
```
For example, the zones in region us-east-2 are `"us-east-2a us-east-2b us-east-2c"`.


#### GCP

To determine available zones in GCP, run the following; available zones within a region can be found in the first column of the output:

```bash
gcloud compute zones list | grep <region>
```

For example, the zones in region us-central1 are `"us-central1-a us-central1-b us-central1c us-central1f"`.


### Determining the AWS batch queue ARN

**Note that if you are using a `miniwdl` engine, you can skip these steps; workflows run via miniwdl will run exclusively in the job queue to which they are submitted.**

1. Visit [the AWS console](https://console.aws.amazon.com/).
2. Navigate to the Batch service.
3. In the lefthand sidebar, select "Compute environments". Note the name of the compute environment with the provisioning model SPOT (if you have deployed a context using spot instances) and the name of the compute environment with provisioning model "EC2" (if you have deployed a context that does not use spot instances).
4. In the lefthand sidebar, select "Job queues".
5. Clicking into an individual queue will show information about the compute environment ("Compute environment order"). Identify the job queue with the Compute environment name that matches the name you identified for the SPOT compute environment; copy the Amazon Resource Name (ARN) for this job queue. This is the value that should be used for the `aws_spot_queue_arn`. Repeat this process to find the ARN for the `aws_on_demand_queue_arn`.

- If `preemptible = true`, only the `aws_spot_queue_arn` is required.
- If `preemptible = false`, only the `aws_on_demand_queue_arn` is required.

# Workflow outputs

## Sample analysis

These files will be output for each sample defined in the cohort.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Array[Array[File]] | bam_stats | Statistics for the set of movie bams for each sample | |
| Array[Array[File]] | read_length_summary | Read length stats for the set of movie bams for each sample | |
| Array[Array[File]] | read_length_quality_summary | Read length quality summaries for the set of movie bams for each sample | |
| Array[Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)]] | aligned_bams | Set of aligned bams for the set of movie bams for each sample | |
| Array[Array[File]] | aligned_bam_mosdepth_summary | Mosdepth summary for the set of aligned bams for each sample | |
| Array[Array[File]] | aligned_bam_mosdepth_region_bed | Mosdepth region bed for the set of aligned bams for each sample | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | small_variant_vcfs | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | small_variant_gvcfs | | |
| Array[File] | small_variant_vcf_stats | | |
| Array[File] | small_variant_roh_bed | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_sv_vcfs | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_phased_small_variant_vcfs | | |
| Array[File] | sample_whatshap_stats_gtfs | | |
| Array[File] | sample_whatshap_stats_tsvs | | |
| Array[File] | sample_whatshap_stats_blocklists | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | merged_haplotagged_bam | | |
| Array[File] | haplotagged_bam_mosdepth_summary | | |
| Array[File] | haplotagged_bam_mosdepth_region_bed | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | trgt_spanning_reads | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | trgt_repeat_vcf | | |
| Array[File] | trgt_dropouts | | |
| Array[Array[File]] | cpg_pileups | | |


## Cohort analysis

These files will be output if the cohort includes more than one sample.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_sv_vcf | | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_phased_joint_called_vcf | | |
| File? | cohort_whatshap_stats_gtfs | | |
| File? | cohort_whatshap_stats_tsvs | | |
| File? | cohort_whatshap_stats_blocklists | | |


## De novo assembly - sample

These files will be output if `cohort.samples[sample]` is set to `true` for any sample.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Array[Array[File]?] | assembly_noseq_gfas | | |
| Array[Array[File]?] | assembly_lowQ_beds | | |
| Array[Array[File]?] | zipped_assembly_fastas | | |
| Array[Array[File]?] | assembly_stats | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)?] | asm_bam | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)?] | htsbox_vcf | | |
| Array[File?] | htsbox_vcf_stats | | |


## De novo assembly - trio

These files will be output if `cohort.de_novo_assembly_trio` is set to `true`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Array[Map[String, String]]? | haplotype_key | | |
| Array[Array[File]]? | trio_assembly_noseq_gfas | | |
| Array[Array[File]]? | trio_assembly_lowQ_beds | | |
| Array[Array[File]]? | trio_zipped_assembly_fastas | | |
| Array[Array[File]]? | trio_assembly_stats | | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)]? | trio_asm_bams | | |


## Tertiary analysis

These files will be output for each run of the workflow. The files that are being annotated will depend on whether the number of samples is equal to or greater than one:
- If the number of samples is equal to one, the files being annotated in this step are the sample small variant VCF and SV VCF.
- If the number of samples is greater than one, the files being annotated in this step are the phased, joint-called VCF and the cohort SV VCF.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | filtered_small_variant_vcf | | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | compound_het_small_variant_vcf | | |
| File | filtered_small_variant_tsv | | |
| File | compound_het_small_variant_tsv | | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | filtered_svpack_vcf | | |
| File | filtered_svpack_tsv | | |
