version 1.0

struct Sample {
  ## Unique identifier for the sample
  String sample_id
  ## Sample sex
  String? sex
  ## Whether the sample is affected by the condition being studied
  Boolean affected

  ## hifi_reads BAMs for the sample
  Array[File] hifi_reads
  ## fail_reads BAMs for the sample
  Array[File]? fail_reads

  ## Optional identifiers for the sample's father, if known
  String? father_id
  ## Optional identifiers for the sample's mother, if known
  String? mother_id
}

struct Family {
  ## Unique identifier for the family
  String family_id
  ## Samples in the family
  Array[Sample] samples
}
