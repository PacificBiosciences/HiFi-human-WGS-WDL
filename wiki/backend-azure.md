# Configuring Cromwell on Azure

Workflows can be run in Azure by setting up [Cromwell on Azure (CoA)](https://github.com/microsoft/CromwellOnAzure). Documentation on deploying and configuring an instance of CoA can be found [here](https://github.com/microsoft/CromwellOnAzure/wiki/Deploy-your-instance-of-Cromwell-on-Azure).

## Requirements

- [Cromwell on Azure](https://github.com/microsoft/CromwellOnAzure) version 3.2+; version 4.0+ is recommended

## Configuring and running the workflow

### Filling out workflow inputs

Fill out any information missing in [the inputs file](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/singleton.azure.inputs.json).

See [the inputs section of the main README](./singleton#inputs) for more information on the structure of the inputs.json file.

### Running via Cromwell on Azure

Instructions for running a workflow from Cromwell on Azure are described in [the Cromwell on Azure documentation](https://github.com/microsoft/CromwellOnAzure/wiki/Running-Workflows).

## Reference data hosted in Azure

To use Azure reference data, add the following line to your `containers-to-mount` file in your Cromwell on Azure installation ([more info here](https://github.com/microsoft/CromwellOnAzure/blob/develop/docs/troubleshooting-guide.md#use-input-data-files-from-an-existing-azure-storage-account-that-my-lab-or-team-is-currently-using)):

`https://datasetpbrarediseases.blob.core.windows.net/dataset?si=public&spr=https&sv=2021-06-08&sr=c&sig=o6OkcqWWlGcGOOr8I8gCA%2BJwlpA%2FYsRz0DMB8CCtCJk%3D`

The [Azure input file template](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/singleton.azure.inputs.json) has paths to the reference files in this blob storage prefilled.
