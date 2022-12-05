# wdl-humanwgs

PacBio's human WGS workflow written in [Workflow Description Language (WDL)](https://openwdl.org/).

For the snakemake version of these workflows, see [here](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake).

Docker images used by these workflows are defined [here](https://github.com/PacificBiosciences/wdl-dockerfiles).


## Inputs

### `Cohort cohort`

A cohort can include one or more samples. Samples need not be related.

- `String cohort_id`: A unique name for the cohort; used to name outputs
- `Array[Sample] samples`: The set of samples for the cohort; see [samples](#samples)
- `File pedigree`: The pedigree for the cohort


#### Samples

- `String sample_id`: A unique name for the sample; used to name outputs
- `Array[IndexData] movie_bams`: The set of movie bams associated with this sample
- `Array[String] phenotypes`: [HPO phenotypes](https://hpo.jax.org/app/) associated with the sample
- `Boolean affected`: The affected status for the sample
- `String? mother_id`: Maternal `sample_id`, if available
- `String? father_id`: Paternal `sample_id`, if available


### `ReferenceData reference`

- `String name`: Reference name; used to name outputs
- `IndexData reference_genome`: Reference genome and index to align reads to
- `Array[String] chromosomes`: Chromosomes to phase during WhatsHap phasing
- `File chromosome_lengths`: File specifying the lengths of each of the reference chromosomes
- `File tandem_repeat_bed`: Tandem repeat locations in the reference genome
- `File trgt_tandem_repeat_bed`: Repeat bed used by TRGT to output spanning reads and a repeat VCF
- `File gnomad_af`: gnomAD allele frequencies; used for annotation
- `File hprc_af`: Allele frequences from the Human Pangenome Reference Consortium
- `File gff`: GFF3 annotation file


### Other inputs

- `File slivar_js`: Additional javascript functions for slivar
- `String deepvariant_version`: Version of deepvariant to use
- `File? deepvariant_model`: Optional alternate deepvariant model file to use
- `Boolean run_de_novo_assembly`: Run the de novo assembly pipeline [false]
- `String container_registry`: Container registry where docker images are hosted


## Workflows

### Main

Calls all subworkflows for a full analysis.

**Workflow**: [workflows/main.wdl](workflows/main.wdl)
**Inputs**: [workflows/inputs.json](workflows/inputs.json)


### [Common](workflows/common)

These are resources that are used across workflows, including structs and common tasks.


### SMRT Cell analysis

Aligns reads to a reference genome and generates statistics on alignment depth, read length, and alignment quality.

**Workflow**: [workflows/smrtcell_analysis/smrtcell_analysis.wdl](workflows/smrtcell_analysis/smrtcell_analysis.wdl)
**Inputs**: [workflows/smrtcell_analysis/inputs.json](workflows/smrtcell_analysis/inputs.json)


### Sample analysis

Calls and phases small and large variants and structural variants.

**Workflow**: [workflows/sample_analysis/sample_analysis.wdl](workflows/sample_analysis/sample_analysis.wdl)
**Inputs**: [workflows/sample_analysis.inputs.json](workflows/sample_analysis/inputs.json)


### Single-sample _de novo_ assembly

Assembles a single genome.

**Workflow**: [workflows/de_novo_assembly/de_novo_assembly.wdl](workflows/de_novo_assembly/de_novo_assembly.wdl)
**Inputs**: [workflows/de_novo_assembly/inputs.json](workflows/de_novo_assembly/inputs.json)


### Cohort analysis

Runs joint genotyping and annotation for a cohort.

**Workflow**: [workflows/cohort_analysis/cohort_analysis.wdl](workflows/cohort_analysis/cohort_analysis.wdl)
