# Configuring Cromwell on GCP

[Cromwell's documentation](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/) on getting started with Google's genomics Pipelines API can be used to set up the resources needed to run the workflow.

## Configuring and running the workflow

### Filling out workflow inputs

Fill out any information missing in [the inputs file](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/gcp/singleton.gcp.inputs.json).

See [the inputs section of the singleton README](./singleton.md#inputs) for more information on the structure of the inputs.json file.

#### Determining available zones

To determine available zones in GCP, run the following; available zones within a region can be found in the first column of the output:

```bash
gcloud compute zones list | grep <region>
```

For example, the zones in region `us-central1` are `"us-central1-a us-central1-b us-central1c us-central1f"`.

## Running the workflow via Google's genomics Pipelines API

[Cromwell's documentation](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/) on getting started with Google's genomics Pipelines API can be used as an example for how to run the workflow.

## Reference data hosted in GCP

GCP reference data is hosted in the `us-west1` region in the bucket `gs://pacbio-wdl`. This bucket is requester-pays, meaning that users will need to [provide a billing project in their Cromwell configuration](https://cromwell.readthedocs.io/en/stable/filesystems/GoogleCloudStorage/) in order to use files located in this bucket.

To avoid egress charges, Cromwell should be set up to spin up compute resources in the same region in which the data is located. If possible, add cohort data to the same region as the reference dataset, or consider mirroring this dataset in the region where your data is located. See [Google's information about data storage and egress charges for more information](https://cloud.google.com/storage/pricing).
