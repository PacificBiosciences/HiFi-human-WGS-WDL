# Configuring Cromwell on Azure

Workflows can be run in Azure by setting up [Cromwell on Azure (CoA)](https://github.com/microsoft/CromwellOnAzure). Documentation on deploying and configuring an instance of CoA can be found [here](https://github.com/microsoft/CromwellOnAzure/wiki/Deploy-your-instance-of-Cromwell-on-Azure).

## Requirements

- [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) version 3.2+; version 4.0+ is recommended

## Configuring and running the workflow

### Filling out workflow inputs

Fill out any information missing in [the inputs file](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/singleton.azure.inputs.json).

See [the inputs section of the main README](./singleton.md#inputs) for more information on the structure of the inputs.json file.

### Running via Cromwell on Azure

Instructions for running a workflow from Cromwell on Azure are described in [the Cromwell on Azure documentation](https://github.com/microsoft/CromwellOnAzure/wiki/Running-Workflows).

## Reference data

Reference data must be uploaded to Azure Blob Storage

[<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17086906.svg" alt="10.5281/zenodo.17086906">](https://zenodo.org/records/17086906)

Reference data is hosted on Zenodo at [10.5281/zenodo.17086906](https://zenodo.org/record/17086906). Download the reference data bundle and extract it to a location on your HPC, then update the input template file with the path to the reference data.

```bash
## download the reference data bundle
wget https://zenodo.org/record/17086906/files/hifi-wdl-resources-v3.1.0.tar

## extract the reference data bundle and rename as dataset
tar -xvf hifi-wdl-resources-v3.1.0.tar
```

Upload the hifi-wdl-resources-v3.1.0 directory to Azure Blob Storage. The reference data can be uploaded to any blob storage container for which your Cromwell on Azure instance has access.  Update the [Azure input file template](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/singleton.azure.inputs.json) with the path to the reference data in blob storage.
