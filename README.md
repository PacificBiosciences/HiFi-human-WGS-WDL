<h1 align="center"><img width="300px" src="https://raw.githubusercontent.com/PacificBiosciences/HiFi-human-WGS-WDL/refs/heads/main/images/logo_wdl_workflows.svg" alt="PacBio WGS Variant Pipeline"/></h1>

<h1 align="center">PacBio WGS Variant Pipeline</h1>

Workflow for analyzing human PacBio whole genome sequencing (WGS) data using [Workflow Description Language (WDL)](https://openwdl.org/).

- Docker images used by this workflow are defined in [the wdl-dockerfiles repo](../../../wdl-dockerfiles). Images are hosted in PacBio's [quay.io repo](https://quay.io/organization/pacbio).

## Workflow

Starting in v2, this repo contains two related workflows. The `singleton` workflow is designed to analyze a single sample, while the `family` workflow is designed to analyze a family of related samples. With the exception of the joint calling tasks in the `family` workflow, both workflows make use of the same tasks, although the input and output structure differ.

The `family` workflow will be best for most use cases. The `singleton` workflow inputs and output structures are relatively flat, which should improve compatibility with platforms like Terra.

Both workflows are designed to analyze PacBio human whole genome sequencing (WGS) data. The workflows are primarily tested and designed to be run on AWS HealthOmics and HPC backends, but are most likely still compatible with Azure (Cromwell on Azure) and GCP (Pipelines API).

**Workflow entrypoint**:

- [workflows/singleton.wdl](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/workflows/singleton.wdl)
- [workflows/family.wdl](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/workflows/family.wdl)

## Setup

This is an actively developed workflow with multiple versioned releases, and we make use of git submodules for common tasks that are shared by multiple workflows. There are two ways to ensure you are using a supported release of the workflow and ensure that the submodules are correctly initialized:

1) Download the release zips directly from a [supported release](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/releases/tag/v4.0.0):

   ```bash
   wget https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/releases/download/v4.0.0/hifi-human-wgs-singleton.zip
   wget https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/releases/download/v4.0.0/hifi-human-wgs-family.zip
   ```

2) Clone the repository and initialize the submodules:

   ```bash
   git clone \
     --depth 1 --branch v4.0.0 \
     https://github.com/PacificBiosciences/HiFi-human-WGS-WDL.git
   ```

## Resource requirements

The most compute-heavy step in the workflow requires 64 cpu cores, and the most memory-heavy step uses 128GiB of RAM. Ensure that the backend environment you're using has enough quota to run the workflow.

On some backends, you may be able to make use of a GPU to accelerate the DeepVariant step.The GPU is not required, but it can significantly speed up the workflow. If you have access to a GPU, you can set the `use_gpu` parameter to `true` in the inputs JSON file.

## Setting up and executing the workflow

1. [Select a backend environment](#selecting-a-backend)
2. [Configure a workflow execution engine in the chosen environment](#configuring-a-workflow-engine-and-container-runtime)
3. [Fill out the inputs JSON file for your cohort](#filling-out-the-inputs-json)
4. [Run the workflow](#running-the-workflow)

### Selecting a backend

The workflow can be run on Azure, AWS, GCP, or HPC. Your choice of backend will largely be determined by the location of your data.

For backend-specific configuration, see the relevant documentation:

- [Azure](./docs/backend-azure.md)
- [AWS](./docs/backend-aws-healthomics.md)
- [GCP](./docs/backend-gcp.md)
- [HPC](./docs/backend-hpc.md)

### Configuring a workflow engine and container runtime

An execution engine is required to run workflows. Three popular engines for running WDL-based workflows are [`sprocket`](https://sprocket.bio), [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/getting_started.html) and [`Cromwell`](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

Because workflow dependencies are containerized, a container runtime is required. This workflow has been tested with [Docker](https://docs.docker.com/get-docker/) and [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/) container runtimes.

See the backend-specific documentation for details on setting up an engine.

### Filling out the inputs JSON

The input to a workflow run is defined in JSON format. Template input files with reference dataset information filled out are available for each backend:

- [HPC singleton entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/hpc/singleton.hpc.inputs.json)
- [HPC family entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/hpc/family.hpc.inputs.json)
- [AWS singleton entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/aws-healthomics/singleton.healthomics.inputs.json)
- [AWS family entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/aws-healthomics/family.healthomics.inputs.json)
- [Azure singleton entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/singleton.azure.inputs.json)
- [Azure family entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/family.azure.inputs.json)
- [GCP singleton entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/gcp/singleton.gcp.inputs.json)
- [GCP family entrypoint](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/gcp/family.gcp.inputs.json)

Using the appropriate inputs template file, fill in the cohort and sample information (see [Workflow Inputs](#workflow-inputs) for more information on the input structure).

### Running the workflow

Run the workflow using the engine and backend that you have configured ([miniwdl](#run-directly-using-miniwdl), [Cromwell](#run-directly-using-cromwell)).

Note that the calls to `miniwdl` and `Cromwell` assume you are accessing the engine directly on the machine on which it has been deployed. Depending on the backend you have configured, you may be able to submit workflows using different methods (e.g. using trigger files in Azure, or using the Amazon Genomics CLI in AWS).

## Workflow inputs

This section describes the inputs required for a run of the workflow. Typically, only the sample-specific sections will be filled out by the user for each run of the workflow. Input templates with reference file locations filled out are provided [for each backend](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends).

Workflow inputs for each entrypoint are described in [singleton](./docs/singleton.md) and [family](./docs/family.md) documentation.

Static reference assets (e.g., reference genome FASTA sequence, or BED files used tools) are packaged within a container that is called at the beginning of the workflow. There is no need to download these assets to run the workflow, but the resource bundle containing the GRCh38 references and other files used in this workflow can be downloaded from Zenodo:

[<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.21517827.svg" alt="10.5281/zenodo.21517827">](https://zenodo.org/records/21517827)

## Tool versions and Docker images

Docker images definitions used by this workflow can be found in [the wdl-dockerfiles repository](../../../wdl-dockerfiles/). Images are hosted in PacBio's [quay.io repo](https://quay.io/organization/pacbio). Docker images used in the workflow are pinned to specific versions by referring to their digests rather than tags.

The Docker image used by a particular step of the workflow can be identified by looking at the `docker` key in the `runtime` block for the given task. Images can be referenced in the following table by looking for the name after the final `/` character and before the `@sha256:...`. For example, the image referred to here is "align_hifiasm":
> ~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f ... b70a8e87

Tool versions and Docker images used in these workflows can be found in the [tools and containers](./docs/tools_containers.md) documentation. For a curated list of the primary scientific tools and links to their own documentation, see [key tools](./docs/tools.md).

---

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
