# GPU support

Workflows can optionally run small variant calling on GPU-enabled nodes. Two inputs control this:

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Boolean | use_gpu | Use GPU for DeepVariant's `call_variants` step. | default = `false` |
| Boolean | use_parabricks_deepvariant | Use NVIDIA Parabricks DeepVariant instead of standard DeepVariant. Only takes effect when `use_gpu` is also `true`. | default = `false` |
| String | gpuType | Type of GPU/Accelerator to use. | Required when `use_gpu` is `true` on a cloud backend; must match backend. |

## Where GPU is used

- **[DeepVariant](deepvariant.md)** — when `use_gpu=true` (and Parabricks is not selected), the `call_variants` step runs on 1 GPU/8 CPU/32GiB instead of 64 CPU/256GiB. `make_examples` and `postprocess_variants` always run on CPU.
- **[Parabricks DeepVariant](parabricks.md)** — when `use_gpu=true` and `use_parabricks_deepvariant=true`, small variant calling is replaced entirely by NVIDIA Parabricks' `pbrun deepvariant`, which runs on 4 GPUs/48 CPU/192GiB and requires GPUs with ≥16GB memory (e.g. V100 or A100). If `use_parabricks_deepvariant=true` but `use_gpu=false`, the workflow falls back to standard CPU DeepVariant.

## GPU Types

| Backend | GPU Type | Notes |
| ------- | -------- | ----- |
| GCP | `nvidia-tesla-k80`, `nvidia-tesla-p100`, `nvidia-tesla-v100`, `nvidia-tesla-p4`, `nvidia-tesla-t4`, `nvidia-tesla-a100`, `nvidia-a100-80gb`, `nvidia-l4`, `nvidia-h100-80gb` | [GPU availability varies by zone.](https://cloud.google.com/compute/docs/gpus/gpu-regions-zones) |
| AWS-HealthOmics | `nvidia-tesla-t4`, `nvidia-tesla-t4-a10g`, `nvidia-tesla-a10g`, `nvidia-t4-a10g-l4`, `nvidia-l4-a10g`, `nvidia-l4`, `nvidia-l40s` | Maps to HealthOmics instance families: `nvidia-tesla-t4`→G4, `nvidia-tesla-t4-a10g`→G4/G5, `nvidia-tesla-a10g`→G5, `nvidia-t4-a10g-l4`→G4/G5/G6, `nvidia-l4-a10g`→G5/G6, `nvidia-l4`→G6, `nvidia-l40s`→G6e. **Recommended: `nvidia-t4-a10g-l4`** — it spans the most instance families (G4/G5/G6) and is what we've tested on us-east-1. |
| Azure | | GPUs are not available on Azure. |
| HPC | | Depends on HPC and miniwdl/Cromwell/sprocket configuration. Reach out to [support@pacb.com](mailto:support@pacb.com?subject=WDL%20Workflows%20-%20GPU%20Support) |
