version 1.0

struct IndexData {
	File data
	File data_index
}

struct DeepVariantModel {
	IndexData model
	File metadata
}

struct RuntimeAttributes {
	# The number of times to retry a task that fails due to preemption
	Int preemptible_tries
	# The number of times to retry a task that fails due a to nonzero return code
	Int max_retries

	String zones
	String queue_arn
	String container_registry
}