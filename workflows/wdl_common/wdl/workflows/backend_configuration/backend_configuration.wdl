version 1.0

# Set runtime attributes across environments depending on the backend in use

import "../../structs.wdl"

workflow backend_configuration {
	input {
		String backend
		String? zones
		String? aws_spot_queue_arn
		String? aws_on_demand_queue_arn
		String container_registry = "quay.io/pacbio"
	}

	if (backend == "GCP") {
		# zones must be defined

		# preemptible_tries applies to failures due to preemption only
		# max_retries applies to failures due to a nonzero rc
		# queue_arn is not used in GCP
		RuntimeAttributes gcp_spot_runtime_attributes = {
			"preemptible_tries": 3,
			"max_retries": 0,
			"zones": select_first([zones]),
			"queue_arn": "",
			"container_registry": container_registry
		}

		RuntimeAttributes gcp_on_demand_runtime_attributes = {
			"preemptible_tries": 0,
			"max_retries": 0,
			"zones": select_first([zones]),
			"queue_arn": "",
			"container_registry": container_registry
		}
	}

	if (backend == "Azure") {
		# Requires Cromwell on Azure v3.2+
		# preemptible_tries >= 1 will be converted to `true`; 0 will be converted to `false`
		# max_retries applies to failures due to preemption or to a nonzero rc
		# zones, queue_arn not used in Azure
		RuntimeAttributes azure_spot_runtime_attributes = {
			"preemptible_tries": 3,
			"max_retries": 3,
			"zones": "",
			"queue_arn": "",
			"container_registry": container_registry
		}

		RuntimeAttributes azure_on_demand_runtime_attributes = {
			"preemptible_tries": 0,
			"max_retries": 0,
			"zones": "",
			"queue_arn": "",
			"container_registry": container_registry
		}
	}

	if (backend == "AWS") {
		# zones must be defined
		# aws_spot_queue_arn must be defined if preemptible is set to true and engine is not miniwdl
		# aws_on_demand_queue_arn must be defined if preemptible is set to false and engine is not miniwdl
		# Using miniwdl engine, the queue ARN of the context the workflow has been submitted to will be used;
		#   the queue_arn runtime attribute will be ignored

		# max_retries applies to failures due to preemption or to a nonzero rc
		# preemptible is not used in AWS
		RuntimeAttributes aws_spot_runtime_attributes = {
			"preemptible_tries": 3,
			"max_retries": 3,
			"zones": select_first([zones]),
			"queue_arn": select_first([aws_spot_queue_arn, ""]),
			"container_registry": container_registry
		}

		RuntimeAttributes aws_on_demand_runtime_attributes = {
			"preemptible_tries": 0,
			"max_retries": 0,
			"zones": select_first([zones]),
			"queue_arn": select_first([aws_on_demand_queue_arn, ""]),
			"container_registry": container_registry
		}
	}

	if (backend == "HPC") {
		# No distinction between preemptible and on-demand in HPC configuration
		RuntimeAttributes hpc_runtime_attributes = {
			"preemptible_tries": 0,
			"max_retries": 3,
			"zones": "",
			"queue_arn": "",
			"container_registry": container_registry
		}
	}

	output {
		RuntimeAttributes spot_runtime_attributes = select_first([
			gcp_spot_runtime_attributes,
			azure_spot_runtime_attributes,
			aws_spot_runtime_attributes,
			hpc_runtime_attributes
		])
		RuntimeAttributes on_demand_runtime_attributes = select_first([
			gcp_on_demand_runtime_attributes,
			azure_on_demand_runtime_attributes,
			aws_on_demand_runtime_attributes,
			hpc_runtime_attributes
		])
	}

	parameter_meta {
		backend: {help: "Backend where the workflow will be executed ['GCP', 'Azure', 'AWS', 'HPC']"}
		zones: {help: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"}
		aws_spot_queue_arn: {help: "Queue ARN for the spot batch queue; required if backend is set to 'AWS'"}
		aws_on_demand_queue_arn: {help: "Queue ARN for the on demand batch queue; required if backend is set to 'AWS'"}
	}
}
