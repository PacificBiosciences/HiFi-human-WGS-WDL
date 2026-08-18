version 1.0

import "../../structs.wdl"

workflow backend_configuration {
  meta {
    description: "Set runtime attributes depending on the backend in use"
    category: "Utility"
    outputs: {
      spot_runtime_attributes: {
        description: "Runtime attributes for spot/preemptible tasks"
      },
      on_demand_runtime_attributes: {
        description: "Runtime attributes for on-demand tasks"
      }
    }
  }

  parameter_meta {
    backend: {
      description: "Backend where the workflow will be executed",
      choices: [
        "GCP",
        "Azure",
        "AWS-HealthOmics",
        "HPC"
      ]
    }
    zones: {
      description: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      description: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      description: "Type of GPU/Accelerator to use"
    }
    container_registry: {
      description: "Container registry to use"
    }
  }

  input {
    String backend
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry
  }

  String default_container_registry = "quay.io/pacbio"

  if (backend == "GCP") {
    # zones must be defined
    # cpuPlatform may be defined

    # preemptible_tries applies to failures due to preemption only
    # max_retries applies to failures due to a nonzero rc
    # gpuCount and gpuType are optional
    # gpuType: ["nvidia-tesla-k80", "nvidia-tesla-p100", "nvidia-tesla-v100", "nvidia-tesla-p4", "nvidia-tesla-t4",
    #           "nvidia-tesla-a100", "nvidia-a100-80gb", "nvidia-l4", "nvidia-h100-80gb"]
    # TODO: Which are compatible with machine type for deepvariant_call_variants?
    RuntimeAttributes gcp_spot_runtime_attributes = object {
      backend: "GCP",
      preemptible_tries: 3,
      max_retries: 3,
      zones: select_first([
        zones
      ]),
      cpuPlatform: select_first([
        cpuPlatform,
        ""
      ]),
      gpuType: select_first([
        gpuType,
        ""
      ]),
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }

    RuntimeAttributes gcp_on_demand_runtime_attributes = object {
      backend: "GCP",
      preemptible_tries: 0,
      max_retries: 0,
      zones: select_first([
        zones
      ]),
      cpuPlatform: select_first([
        cpuPlatform,
        ""
      ]),
      gpuType: select_first([
        gpuType,
        ""
      ]),
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }
  }

  if (backend == "Azure") {
    # Requires Cromwell on Azure v3.2+
    # preemptible_tries >= 1 will be converted to `true`; 0 will be converted to `false`
    # max_retries applies to failures due to preemption or to a nonzero rc
    # GPUs are not available in Azure
    RuntimeAttributes azure_spot_runtime_attributes = object {
      backend: "Azure",
      preemptible_tries: 1,
      max_retries: 3,
      zones: "",
      cpuPlatform: "",
      gpuType: "",
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }

    RuntimeAttributes azure_on_demand_runtime_attributes = object {
      backend: "Azure",
      preemptible_tries: 0,
      max_retries: 0,
      zones: "",
      cpuPlatform: "",
      gpuType: "",
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }
  }

  if (backend == "AWS-HealthOmics") {
    # preemptible retries AWS service errors
    # maxRetries retries OOM errors with 2x memory, requires GNU findutils 4.2.3+
    # https://docs.aws.amazon.com/omics/latest/dev/workflow-languages-wdl.html#workflow-wdl-directives

    # gpuCount and gpuType are optional
    # tested on us-east-1 with "nvidia-t4-a10g-l4"
    # gpuType options: ["nvidia-tesla-t4", "nvidia-tesla-t4-a10g", "nvidia-tesla-a10g", "nvidia-t4-a10g-l4", "nvidia-l4-a10g", "nvidia-l4", "nvidia-l40s"]
    # Accelerator spec	->	Healthomics instance types
    # "nvidia-tesla-t4"	->	G4
    # "nvidia-tesla-t4-a10g"	->	G4 and G5
    # "nvidia-tesla-a10g"	->	G5
    # "nvidia-t4-a10g-l4"	->	G4, G5, and G6
    # "nvidia-l4-a10g"	->	G5 and G6
    # "nvidia-l4"	->	G6
    # "nvidia-l40s"	->	G6e

    # AWS HealthOmics must use containers hosted on ECR.
    # Our recommended approach for HealthOmics is:
    # 1. Configure ECR pull through cache rules to pull from quay.io (and docker-hub if desired) https://docs.aws.amazon.com/omics/latest/dev/workflows-ecr.html#ecr-pull-through-configure
    # 2. Configure registry mappings to remap quay.io and docker-hub URLs to your ECR https://docs.aws.amazon.com/omics/latest/dev/workflows-ecr.html#ecr-pull-through-registry-mapping
    # 3. Manually mirror NVIDIA Clara containers to your ECR
    # 4. Configure image mappings to remap nvidia clara urls to your ECR https://docs.aws.amazon.com/omics/latest/dev/workflows-ecr.html#ecr-pull-through-mapping-format
    # Check the workflow documentation for more detailed instructions.
    # If you don't plan to configure ECR and mappings, mirror the containers you need to your own ECR and set container_registry to your ECR URL prefix
    RuntimeAttributes aws_healthomics_on_demand_runtime_attributes = object {
      backend: "AWS-HealthOmics",
      preemptible_tries: 2,
      max_retries: 1,
      zones: "",
      cpuPlatform: "",
      gpuType: select_first([
        gpuType,
        ""
      ]),
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }
  }

  if (backend == "HPC") {
    # No distinction between preemptible and on-demand in HPC configuration
    RuntimeAttributes hpc_runtime_attributes = object {
      backend: "HPC",
      preemptible_tries: 0,
      max_retries: 3,
      zones: "",
      cpuPlatform: "",
      gpuType: select_first([
        gpuType,
        ""
      ]),
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }
  }

  output {
    RuntimeAttributes spot_runtime_attributes = select_first([
      gcp_spot_runtime_attributes,
      azure_spot_runtime_attributes,
      aws_healthomics_on_demand_runtime_attributes,
      hpc_runtime_attributes
    ])
    RuntimeAttributes on_demand_runtime_attributes = select_first([
      gcp_on_demand_runtime_attributes,
      azure_on_demand_runtime_attributes,
      aws_healthomics_on_demand_runtime_attributes,
      hpc_runtime_attributes
    ])
  }
}

