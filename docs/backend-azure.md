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

## Reference data hosted in Azure

Azure reference data can be used without any modifications to your Cromwell on Azure by using the templates provided in this repository.

The [Azure input file template](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/azure/singleton.azure.inputs.json) has paths to the reference files in this blob storage prefilled.
