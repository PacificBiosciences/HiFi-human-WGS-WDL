# Parabricks DeepVariant subworkflow

If both `parabricks` and `gpu` are set to `true`, Parabricks DeepVariant will be used instead of standard DeepVariant.  So far, we have only tested this on HPC and AWS backends with NVIDIA GPUs.  The Parabricks DeepVariant implementation requires 4 GPUs with ≥16GB memory (e.g. NVIDIA V100 or A100).  The `parabricks_deepvariant` task will run on 4GPU/48CPU/192GiB.
