# Configuring Cromwell on GCP

[Cromwell's documentation](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/) on getting started with Google's genomics Pipelines API can be used to set up the resources needed to run the workflow.

# Configuring and running the workflow

## Filling out workflow inputs

Fill out any information missing in [the inputs file](inputs.gcp.json).

See [the inputs section of the main README](../../README.md#workflow-inputs) for more information on the structure of the inputs.json file.

### Determining available zones

To determine available zones in GCP, run the following; available zones within a region can be found in the first column of the output:

```bash
gcloud compute zones list | grep <region>
```

For example, the zones in region us-central1 are `"us-central1-a us-central1-b us-central1c us-central1f"`.

## Reference data hosted in GCP

\<TODO>
