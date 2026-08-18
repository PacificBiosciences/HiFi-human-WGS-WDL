version 1.0

## Runtime attributes to be used by tasks and subworkflows.
struct RuntimeAttributes {
  ## Backend where the workflow will be executed
  String backend

  ## The number of times to retry a task that fails due to preemption
  Int preemptible_tries
  ## The number of times to retry a task that fails due a to nonzero return code
  Int max_retries

  ## Zones where compute will take place; required if backend is set to 'GCP'
  String zones
  ## Optional minimum CPU platform to use for tasks on GCP
  String cpuPlatform

  ## Optional type of GPU/Accelerator to use
  String gpuType

  ## Container registry to use, default is "quay.io/pacbio"
  String container_registry
}

