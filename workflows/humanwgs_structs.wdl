version 1.0

struct Sample {
  String sample_id

  String? sex
  Boolean affected

  Array[File] hifi_reads

  String? father_id
  String? mother_id
}

struct Family {
  String family_id
  Array[Sample] samples
}
