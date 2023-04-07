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

```yaml
data:
  - location: s3://dnastack-resources
    readOnly: true
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
| Array[String] | phenotypes | [Human Phenotype Ontology (HPO) phenotypes](https://hpo.jax.org/app/) associated with the cohort | |
| Boolean | run_de_novo_assembly_trio | Run trio binned _de novo_ assembly. | Cohort must contain at least one valid trio (child and both parents present in the cohort) |

### [Sample](workflows/humanwgs_structs.wdl)

Sample information for each sample in the workflow run.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | sample_id | A unique name for the sample; used to name outputs | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | movie_bams | The set of unaligned movie BAMs associated with this sample | |
| String | sex | Sample sex | \["MALE", "FEMALE"\] |
| Boolean | affected | Is this sample affected by the phenotype? | \[true, false\] |
| String? | father_id | Paternal `sample_id` | |
| String? | mother_id | Maternal `sample_id` | |
| Boolean | run_de_novo_assembly | If true, run single-sample _de novo_ assembly for this sample | \[true, false\] |

## [ReferenceData](workflows/humanwgs_structs.wdl)

Files associated with the reference genome.

These files are hosted publicly in each of the cloud backends; see `workflows/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String | name | Reference name; used to name outputs (e.g., "GRCh38") | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl) | fasta | Reference genome and index | |
| Array[String] | chromosomes | Chromosomes to phase, typically `chr{1..22} chr{X,Y}` | |
| File | chromosome_lengths | Reference chromosome lengths | |
| File | tandem_repeat_bed | Tandem repeat locations used by [pbsv](https://github.com/PacificBiosciences/pbsv) to normalize SV representation | |
| File | trgt_tandem_repeat_bed | Tandem repeat sites to be genotyped by [TRGT](https://github.com/PacificBiosciences/trgt) | |
| File | gnomad_af | [gnomAD](https://gnomad.broadinstitute.org/) v3.1 allele frequences in [`slivar gnotate`](https://github.com/brentp/slivar/wiki/gnotate) format | |
| File | hprc_af | Allele frequences in ~100 [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/) samples in [`slivar gnotate`](https://github.com/brentp/slivar/wiki/gnotate) format | |
| File | gff | [Ensembl](https://useast.ensembl.org/index.html) GFF3 reference annotation | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | population_vcfs | An array of structural variant population VCFs | |

## [SlivarData](workflows/humanwgs_structs.wdl)

Files associated with `slivar` annotation.

These files are hosted publicly in each of the cloud backends; see `workflows/inputs.${backend}.json`.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| File | slivar_js | Additional javascript functions for slivar | |
| File | hpo_terms | [HPO](https://hpo.jax.org/app/) annotation lookups | |
| File | hpo_dag | [HPO](https://hpo.jax.org/app/) annotation lookups | |
| File | hpo_annotations | [HPO](https://hpo.jax.org/app/) annotation lookups | |
| File | ensembl_to_hgnc | Ensembl to HGNC gene mapping | |
| File | lof_lookup | Loss-of-function scores per gene | |
| File | clinvar_lookup | ClinVar annotations per gene | |

## Other inputs

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| String? | deepvariant_version | Version of deepvariant to use \["1.5.0"\] | |
| [DeepVariantModel](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | deepvariant_model | Optonal alternate DeepVariant model file to use | |
| Boolean? | run_tertiary_analysis | Run the optional tertiary analysis steps \[true\] | |
| String | backend | Backend where the workflow will be executed | \["Azure", "AWS", "GCP"\] |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'. | [Determining available zones in AWS and GCP](#determining-available-zones-in-aws-and-gcp). |
| String? | aws_spot_queue_arn | Queue ARN for the spot batch queue; required if backend is set to 'AWS' and `preemptible` is set to `true` | [Determining the AWS queue ARN](#determining-the-aws-batch-queue-arn) |
| String? | aws_on_demand_queue_arn | Queue ARN for the on demand batch queue; required if backend is set to 'AWS' and `preemptible` is set to `false` | [Determining the AWS queue ARN](#determining-the-aws-batch-queue-arn) |
| Boolean | preemptible | If set to `true`, run tasks preemptibly where possible. On-demand VMs will be used only for tasks that run for >24 hours if the backend is set to GCP. If set to `false`, on-demand VMs will be used for every task. | \[true, false\] |

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
| Array[Array[File]] | bam_stats | TSV of length and quality for each read (per input BAM) | |
| Array[Array[File]] | read_length_summary | For each input BAM, read length distribution (per input BAM) | |
| Array[Array[File]] | read_quality_summary | For each input BAM, read quality distribution (per input BAM) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | small_variant_gvcfs | Small variants (SNPs and INDELs < 50bp) gVCFs called by [DeepVariant](https://github.com/google/deepvariant) (with index) | |
| Array[File] | small_variant_vcf_stats | [`bcftools stats`](https://samtools.github.io/bcftools/bcftools.html#stats) summary statistics for small variants | |
| Array[File] | small_variant_roh_bed | Regions of homozygosity determiend by [`bcftools roh`](https://samtools.github.io/bcftools/howtos/roh-calling.html) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_sv_vcfs | Structural variants called by [pbsv](https://github.com/PacificBiosciences/pbsv) (with index) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | sample_phased_small_variant_vcfs | Small variants called by [DeepVariant](https://github.com/google/deepvariant) and phased by [WhatsHap](https://whatshap.readthedocs.io/en/latest/) (with index) | |
| Array[File] | sample_whatshap_stats_tsvs | Phase block statistics written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats) | |
| Array[File] | sample_whatshap_stats_gtfs | Phase block GTF written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats)  | |
| Array[File] | sample_whatshap_stats_blocklists | Haplotype block list written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | merged_haplotagged_bam | Aligned (by [pbmm2](https://github.com/PacificBiosciences/pbmm2)), haplotagged (by [`whatshap haplotag`](https://whatshap.readthedocs.io/en/latest/guide.html#visualizing-phasing-results)) reads (with index) | |
| Array[File] | haplotagged_bam_mosdepth_summary | [mosdepth](https://github.com/brentp/mosdepth) summary of median depths per chromosome | |
| Array[File] | haplotagged_bam_mosdepth_region_bed | [mosdepth](https://github.com/brentp/mosdepth) BED of median coverage depth per 500 bp window | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | trgt_repeat_vcf | Tandem repeat genotypes from [TRGT](https://github.com/PacificBiosciences/trgt/blob/main/docs/vcf_files.md) (with index) | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)] | trgt_spanning_reads | Fragments of HiFi reads spanning loci genotyped by [TRGT](https://github.com/PacificBiosciences/trgt/blob/main/docs/tutorial.md) (with index) | |
| Array[File] | trgt_dropouts | Regions with insufficient coverage to genotype | |
| Array[Array[File]] | cpg_pileups | 5mCpG site methylation probability pileups generated by [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools#output-files) | |
| Array[File] | paraphase_output | Output generated by [Paraphase](https://github.com/PacificBiosciences/paraphase) | |
| Array[IndexData] | paraphase_realigned_bam | Realigned BAM for selected medically relevant genes in segmental duplications (with index), generated by [Paraphase](https://github.com/PacificBiosciences/paraphase) | |
| Array[Array[File]] | paraphase_vcfs | Phased Variant calls for selected medically relevant genes in segmental duplications, generated by [Paraphase](https://github.com/PacificBiosciences/paraphase) | |

## De novo assembly - sample

These files will be output if `cohort.samples[sample]` is set to `true` for any sample.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Array[Array[File]?] | zipped_assembly_fastas | [_De novo_ dual assembly](http://lh3.github.io/2021/10/10/introducing-dual-assembly) generated by [hifiasm](https://github.com/chhylp123/hifiasm) | |
| Array[Array[File]?] | assembly_noseq_gfas | Assembly graphs in [GFA format](https://github.com/chhylp123/hifiasm/blob/master/docs/source/interpreting-output.rst). | |
| Array[Array[File]?] | assembly_lowQ_beds | Coordinates of low quality regions in BED format. | |
| Array[Array[File]?] | assembly_stats | Assembly size and NG50 stats generated by [calN50](https://github.com/lh3/calN50). | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)?] | asm_bam | [minimap2](https://github.com/lh3/minimap2) alignment of assembly to reference. | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)?] | htsbox_vcf | Naive pileup variant calling of assembly against reference with [`htsbox`](https://github.com/lh3/htsbox) | |
| Array[File?] | htsbox_vcf_stats | [`bcftools stats`](https://samtools.github.io/bcftools/bcftools.html#stats) summary statistics for `htsbox` variant calls | |

## Cohort analysis

These files will be output if the cohort includes more than one sample.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_sv_vcf | Structural variants joint-called by [pbsv](https://github.com/PacificBiosciences/pbsv) (with index) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | cohort_phased_joint_called_vcf | Small variants called by [DeepVariant](https://github.com/google/deepvariant), joint-called by [GLnexus](https://github.com/dnanexus-rnd/GLnexus), and phased by [WhatsHap](https://whatshap.readthedocs.io/en/latest/) (with index) | |
| File? | cohort_whatshap_stats_tsvs | Phase block statistics written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats)  | |
| File? | cohort_whatshap_stats_gtfs | Phase block GTF written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats) | |
| File? | cohort_whatshap_stats_blocklists | Haplotype block list written by [`whatshap stats`](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-stats) | |

## De novo assembly - trio

These files will be output if `cohort.de_novo_assembly_trio` is set to `true` and there is at least one parent-parent-kid trio in the cohort.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| Array[Array[File]]? | trio_zipped_assembly_fastas | [Haplotype-resolved _de novo_ assembly](http://lh3.github.io/2021/10/10/introducing-dual-assembly) of the trio kid generated by [hifiasm](https://github.com/chhylp123/hifiasm) with [trio binning](https://github.com/chhylp123/hifiasm#trio-binning) | |
| Array[Array[File]]? | trio_assembly_noseq_gfas | Assembly graphs in [GFA format](https://github.com/chhylp123/hifiasm/blob/master/docs/source/interpreting-output.rst). | |
| Array[Array[File]]? | trio_assembly_lowQ_beds | Coordinates of low quality regions in BED format. | |
| Array[Array[File]]? | trio_assembly_stats | Assembly size and NG50 stats generated by [calN50](https://github.com/lh3/calN50). | |
| Array[[IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)]? | trio_asm_bams | [minimap2](https://github.com/lh3/minimap2) alignment of assembly to reference. | |
| Array[Map[String, String]]? | haplotype_key | Indication of which haplotype (`hap1`/`hap2`) corresponds to which parent. | |

## Tertiary analysis

These files will be output for each run of the workflow if `run_tertiary_analysis` is set to `true` (this is the default). The files that are being annotated will depend on whether the number of samples is equal to or greater than one:
- If the number of samples is equal to one, the files being annotated in this step are the sample small variant VCF and SV VCF.
- If the number of samples is greater than one, the files being annotated in this step are the phased, joint-called small variant VCF and the cohort SV VCF.

| Type | Name | Description | Notes |
| :- | :- | :- | :- |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | filtered_small_variant_vcf | Small variant calls that are filtered based on population frequency and annotated with cohort information, population frequency, gene, functional impact, etc., using [slivar](https://github.com/brentp/slivar) and [`bcftools csq`](https://samtools.github.io/bcftools/howtos/csq-calling.html) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | compound_het_small_variant_vcf | Compound heterozygotes annotated with cohort information, population frequency, gene, functional impact, etc., using [slivar](https://github.com/brentp/slivar) and [`bcftools csq`](https://samtools.github.io/bcftools/howtos/csq-calling.html) | |
| File? | filtered_small_variant_tsv | Filtered VCFs are reformatted as a human-readable TSV by [`slivar tsv`](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF) | |
| File? | compound_het_small_variant_tsv | Filtered VCFs are reformatted as a human-readable TSV by [`slivar tsv`](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF) | |
| [IndexData](https://github.com/PacificBiosciences/wdl-common/blob/main/wdl/structs.wdl)? | filtered_svpack_vcf | Structural variant calls that are filtered based on population frequency and annotated with cohort information, population frequency, gene, functional impact, etc., using [svpack](https://github.com/PacificBiosciences/svpack) | |
| File? | filtered_svpack_tsv | Filtered VCFs are reformatted as a human-readable TSV by [`slivar tsv`](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF) | |
