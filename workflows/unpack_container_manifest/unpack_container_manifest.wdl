version 1.0

import "../wdl-common/wdl/structs.wdl"

task unpack_container_manifest {
  meta {
    description: "Unpack reference assets from a container manifest."
    category: "Subworkflow"
    outputs: {
      ref_name: "Short name for reference",
      ref_fasta: "Reference genome FASTA",
      ref_index: "Reference genome FASTA index",
      max_norm_female_chrY_depth: "Maximum expected normalized chrY depth for samples without chrY",
      trgt_tandem_repeat_bed: "Tandem Repeat catalog (BED) for TRGT genotyping",
      sawfish_exclude_bed: "Regions to be excluded for Sawfish CNV calls in gzipped BED format",
      sawfish_exclude_bed_index: "Regions to be excluded for Sawfish CNV calls in gzipped BED format index file",
      sawfish_expected_bed_male: "Expected allosome copy number BED for XY samples",
      sawfish_expected_bed_female: "Expected allosome copy number BED for XX samples",
      methbat_region_tsv: "Regions for MethBat methylation profiling in tab-separated format",
      paraphase_genome_build: "Genome reference build parameter for Paraphase",
      run_starphase: "Whether to run StarPhase task",
      manifest_json: "The manifest JSON file describing this container's static workflow inputs."
    }
  }

  parameter_meta {
    unpack_image: {
      description: "Container image containing static resources and unpack script."
    }
    runtime_attributes: {
      description: "Default Runtime Attributes"
    }
  }

  input {
    String unpack_image
    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb = 4
  Int disk_size = 20
  Int boot_disk_size = 20

  # Path to unpack_container.py is stable and baked into the container image.
  String unpack_script = "/opt/scripts/unpack_container.py"
  String container_manifest = "/opt/manifests/manifest.json"

  command <<<
    set -euo pipefail

    python3 "~{unpack_script}" \
      --manifest "~{container_manifest}" \
      --output-dir "."
  >>>

  output {
    String ref_name = read_string("ref_name")
    File ref_fasta = glob("ref_fasta/*")[0]
    File ref_index = glob("ref_index/*")[0]
    Float max_norm_female_chrY_depth = read_float("max_norm_female_chrY_depth")
    File trgt_tandem_repeat_bed = glob("trgt_tandem_repeat_bed/*")[0]
    File sawfish_exclude_bed = glob("sawfish_exclude_bed/*")[0]
    File sawfish_exclude_bed_index = glob("sawfish_exclude_bed_index/*")[0]
    File sawfish_expected_bed_male = glob("sawfish_expected_bed_male/*")[0]
    File sawfish_expected_bed_female = glob("sawfish_expected_bed_female/*")[0]
    File methbat_region_tsv = glob("methbat_region_tsv/*")[0]
    String paraphase_genome_build = read_string("paraphase_genome_build")
    Boolean run_starphase = read_boolean("run_starphase")
    File manifest_json = "manifest.json"
  }

  runtime {
    docker: "~{unpack_image}"
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: boot_disk_size
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
    cacheable: true  # !UnknownRuntimeKey
  }
}

