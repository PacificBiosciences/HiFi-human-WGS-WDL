version 1.0

struct Sample {
  ## Unique identifier for the sample
  String sample_id

  ## hifi_reads BAMs for the sample
  Array[File] hifi_reads
  ## fail_reads BAMs for the sample
  Array[File]? fail_reads
}

struct Family {
  ## Unique identifier for the family
  String family_id
  ## Samples in the family
  Array[Sample] samples
}

