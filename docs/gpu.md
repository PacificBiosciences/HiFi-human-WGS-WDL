# GPU support

Starting in workflow version 2.0.0, we have added support for running workflows on GPU-enabled nodes. The first task to take advantage of this is the [`deepvariant_call_variants` task](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/workflows/wdl-common/wdl/workflows/deepvariant/deepvariant.wdl) in the DeepVariant workflow, which can use 1 GPU. To run the DeepVariant workflow on a GPU-enabled node, you will need to provide some additional configuration in your inputs JSON file.

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Boolean | gpu | Use GPUs. | default = `false` |
| String | gpuType | Type of GPU/Accelerator to use. | This will depend on your backend configuration. |

## GPU Types

| Backend | GPU Type | Notes |
| ------- | -------- | ----- |
| AWS-HealthOmics | `["nvidia-tesla-a10g", "nvidia-tesla-t4", "nvidia-tesla-t4-a10g"]` | [GPU availability varies by zone.](https://aws.amazon.com/ec2/instance-types) |
| Azure |  | GPU support not yet implemented, but monitoring microsoft/ga4gh-tes#717. |
| GCP | `["nvidia-tesla-t4", "nvidia-tesla-v100"]` | [GPU availability varies by zone.](https://cloud.google.com/compute/docs/gpus/gpu-regions-zones) |
| HPC |  | This will depend on HPC and miniwdl or Cromwell configuration. Reach out to [support@pacb.com](mailto:support@pacb.com?subject=WDL%20Workflows%20-%20GPU%20Support) |
