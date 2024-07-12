<h1 align="center"><img width="300px" src="https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/raw/main/images/logo_wdl_workflows.svg"/></h1>

<h1 align="center">PacBio WGS Variant Pipeline</h1>

Workflow for analyzing human PacBio whole genome sequencing (WGS) data using [Workflow Description Language (WDL)](https://openwdl.org/).

- Docker images used by this workflow are defined in [the wdl-dockerfiles repo](https://github.com/PacificBiosciences/wdl-dockerfiles). Images are hosted in PacBio's [quay.io](https://quay.io/organization/pacbio).
- Common tasks that may be reused within or between workflows are defined in [the wdl-common repo](https://github.com/PacificBiosciences/wdl-common).

## Workflow

Starting in v2, this repo contains two related workflows.  The `singleton` workflow is designed to analyze a single sample, while the `family` workflow is designed to analyze a family of related samples.  With the exception of the joint calling tasks in the `family` workflow, both workflows make use of the same tasks, although the input and output structure differ.

The `singleton` workflow will be best for most use cases, including analysis of rare disease probands without data from relatives, or non-rare disease applications.  The input and output structures are relatively flat, which should improve compatibility with platforms like Terra.

The `family` workflow adds joint-calling for small variants and structural variants and is best suited for rare disease applications where data is available from multiple samples in a family.  The input and output structures are more complex and may require some additional work to integrate with other platforms.

Both workflows are designed to analyze human PacBio whole genome sequencing (WGS) data.  The workflows are designed to be run on Azure, AWS HealthOmics, GCP, or HPC backends.

**Workflow entrypoint**:

- [workflows/singleton.wdl](./workflows/singleton.wdl)
- [workflows/family.wdl](./workflows/family.wdl)

## Setup

This is an actively developed workflow with multiple versioned releases, and we make use of git submodules for common tasks that are shared by multiple workflows. There are two ways to ensure you are using a supported release of the workflow and ensure that the submodules are correctly initialized:

1) Download the release zips directly from a [supported release](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/releases/tag/v2.0.0-rc1):

  ```bash
  wget https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/releases/download/v2.0.0-rc1/hifi-human-wgs-singleton.v2.0.0-rc1.zip
  wget https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/releases/download/v2.0.0-rc1/hifi-human-wgs-family.v2.0.0-rc1.zip
  ```

2) Clone the repository and initialize the submodules:

  ```bash
  git clone \
    --depth 1 --branch v2.0.0-rc1 \
    --recursive \
    https://github.com/PacificBiosciences/HiFi-human-WGS-WDL.git
  ```

## Resource requirements

The most resource-heavy step in the workflow requires 64 cores and 256 GB of RAM. Ensure that the backend environment you're using has enough quota to run the workflow.

On some backends, you may be able to make use of a GPU to accelerate the DeepVariant step.  The GPU is not required, but it can significantly speed up the workflow.  If you have access to a GPU, you can set the `gpu` parameter to `true` in the inputs JSON file.

## Reference datasets and associated workflow files

Reference datasets are hosted publicly for use in the pipeline. For data locations, see the [backend-specific documentation](./docs/backends/) and template inputs files for each backend with paths to publicly hosted reference files filled out.

## Setting up and executing the workflow

1. [Select a backend environment](#selecting-a-backend)
2. [Configure a workflow execution engine in the chosen environment](#configuring-a-workflow-engine)
3. [Fill out the inputs JSON file for your cohort](#filling-out-the-inputs-json)
4. [Run the workflow](#running-the-workflow-1)

### Selecting a backend

The workflow can be run on Azure, AWS, GCP, or HPC. Your choice of backend will largely be determined by the location of your data.

For backend-specific configuration, see the relevant documentation:

- [Azure](./docs/backends/azure)
- [AWS](./docs/backends/aws-healthomics)
- [GCP](./docs/backends/gcp)
- [HPC](./docs/backends/hpc)

### Configuring a workflow engine and container runtime

An execution engine is required to run workflows. Two popular engines for running WDL-based workflows are [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/getting_started.html) and [`Cromwell`](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

Because workflow dependencies are containerized, a container runtime is required. This workflow has been tested with [Docker](https://docs.docker.com/get-docker/) and [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/) container runtimes.

See the [backend-specific documentation](backends) for details on setting up an engine.

| Engine | Azure | AWS | GCP | HPC |
| :- | :- | :- | :- | :- |
| [**miniwdl**](https://github.com/chanzuckerberg/miniwdl#scaling-up) | _Unsupported_ | Supported via [AWS HealthOmics](https://aws.amazon.com/healthomics/) | _Unsupported_ | (SLURM only) Supported via the [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm) plugin |
| [**Cromwell**](https://cromwell.readthedocs.io/en/stable/backends/Backends/) | Supported via [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) | _Unsupported_ | Supported via Google's [Pipelines API](https://cromwell.readthedocs.io/en/stable/backends/Google/) | Supported - [Configuration varies depending on HPC infrastructure](https://cromwell.readthedocs.io/en/stable/tutorials/HPCIntro/) |

### Filling out the inputs JSON

The input to a workflow run is defined in JSON format. Template input files with reference dataset information filled out are available for each backend:

- [Azure](./docs/backends/azure/singleton.azure.inputs.json)
- [AWS](./docs/backends/aws-healthomics/singleton.healthomics.inputs.json)
- [GCP](./docs/backends/gcp/singleton.gcp.inputs.json)
- [HPC](./docs/backends/hpc/singleton.hpc.inputs.json)

Using the appropriate inputs template file, fill in the cohort and sample information (see [Workflow Inputs](#workflow-inputs) for more information on the input structure).

If using an HPC backend, you will need to download the reference bundle and replace the `<local_path_prefix>` in the input template file with the local path to the reference datasets on your HPC.  If using Amazon HealthOmics, you will need to download the reference bundle, upload it to your S3 bucket, and adjust paths accordingly.

### Running the workflow

Run the workflow using the engine and backend that you have configured ([miniwdl](#run-directly-using-miniwdl), [Cromwell](#run-directly-using-cromwell)).

Note that the calls to `miniwdl` and `Cromwell` assume you are accessing the engine directly on the machine on which it has been deployed. Depending on the backend you have configured, you may be able to submit workflows using different methods (e.g. using trigger files in Azure, or using the Amazon Genomics CLI in AWS).

#### Run directly using miniwdl

`miniwdl run workflows/singleton.wdl -i <input_file_path.json>`

#### Run directly using Cromwell

`java -jar <cromwell_jar_path> run workflows/singleton.wdl -i <input_file_path.json>`

If Cromwell is running in server mode, the workflow can be submitted using cURL. Fill in the values of CROMWELL_URL and INPUTS_JSON below, then from the root of the repository, run:

## Workflow inputs

This section describes the inputs required for a run of the workflow. Typically, only the `humanwgs.cohort` and potentially run/backend-specific sections will be filled out by the user for each run of the workflow. Input templates with reference file locations filled out are provided [for each backend](./docs/backends).

Workflow inputs for each entrypoint are described in [singleton](./docs/singleton.md) and [family](./docs/family.md) documentation.

At a high level, we have two types of inputs files:

- _maps_ are TSV files describing inputs that will be used for every execution of the workflow, like reference genome FASTA files and genome interval BED files.
- _inputs.json_ files are JSON files that describe the samples to be analyzed and the paths to the input files for each sample.

The resource bundle containing the GRCh38 reference and other files used in this workflow can be downloaded from Zenodo:

![https://doi.org/10.5281/zenodo.12729255](https://zenodo.org/badge/DOI/10.5281/zenodo.12729255.svg)

# Tool versions and Docker images

Docker images definitions used by this workflow can be found in [the wdl-dockerfiles repository](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/987efde4d614a292fbfe9f3cf146b63005ad6a8a). Images are hosted in PacBio's [quay.io](https://quay.io/organization/pacbio). Docker images used in the workflow are pegged to specific versions by referring to their digests rather than tags.

The Docker image used by a particular step of the workflow can be identified by looking at the `docker` key in the `runtime` block for the given task. Images can be referenced in the following table by looking for the name after the final `/` character and before the `@sha256:...`. For example, the image referred to here is "align_hifiasm":
> ~{runtime_attributes.container_registry}/**align_hifiasm**@sha256:3968cb<...>b01f80fe

Tool versions and Docker images used in these workflows can be found in the [tools_containers](./docs/tools_containers.md) documentation.

---

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
