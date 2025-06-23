version 1.1

struct Sample {
  String sample_id

  String? sex
  Boolean affected

  #Array[File] hifi_reads
  # parental inputs (formerly “normal”); empty if none provided
  Array[File] parental_bams
  # wild-type inputs (formerly “tumor”); must be non-empty
  Array[File] wild_type_bams


  String? father_id
  String? mother_id
}

struct Family {
  String family_id
  Array[Sample] samples
}
