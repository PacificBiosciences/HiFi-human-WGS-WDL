version 1.0

## The `Sample` struct contains sample specific data and metadata.
struct Sample {
  ## Unique identifier for the sample
  String sample_id

  ## Array of paths to hifi_reads in unaligned BAM format
  Array[File] hifi_reads
  ## Array of paths to fail_reads in unaligned BAM format
  Array[File]? fail_reads
}

## The `Family` struct contains the samples for the family.
struct Family {
  ## Unique identifier for the family
  String family_id
  ## Samples in the family
  Array[Sample] samples
}

