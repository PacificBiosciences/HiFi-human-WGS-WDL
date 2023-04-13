## Reference data hosted in Azure

To use Azure reference data, add the following line to your `containers-to-mount` file in your Cromwell on Azure installation ([more info here](https://github.com/microsoft/CromwellOnAzure/blob/develop/docs/troubleshooting-guide.md#use-input-data-files-from-an-existing-azure-storage-account-that-my-lab-or-team-is-currently-using)):

`https://datasetpbrarediseases.blob.core.windows.net/dataset?si=public&spr=https&sv=2021-06-08&sr=c&sig=o6OkcqWWlGcGOOr8I8gCA%2BJwlpA%2FYsRz0DMB8CCtCJk%3D`

The [Azure input file template](inputs.azure.json) has paths to the reference files in this blob storage prefilled.
