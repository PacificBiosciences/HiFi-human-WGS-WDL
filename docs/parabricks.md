# Parabricks

## Parabricks DeepVariant subworkflow

If both `use_parabricks_deepvariant` and `use_gpu` are set to `true`, Parabricks
DeepVariant will be used instead of standard DeepVariant. We have tested this on
HPC and AWS backends with NVIDIA GPUs. The Parabricks DeepVariant implementation
requires 4 GPUs with ≥16GB memory (e.g. NVIDIA V100 or A100). The
`parabricks_deepvariant` task will run on 4GPU/48CPU/192GiB.
