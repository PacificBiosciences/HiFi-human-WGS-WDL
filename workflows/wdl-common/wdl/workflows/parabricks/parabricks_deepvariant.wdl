version 1.0

import "../../structs.wdl"

workflow parabricks_deepvariant {
  meta {
    description: "Call variants from aligned HiFi reads using Parabricks DeepVariant"
    outputs: {
      vcf: {
        description: "Small variant VCF"
      },
      vcf_index: {
        description: "Index for small variant VCF"
      },
      gvcf: {
        description: "Small variant GVCF"
      },
      gvcf_index: {
        description: "Index for small variant GVCF"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Aligned BAM index"
    }
    regions_bed: {
      description: "Regions BED"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    ref_name: {
      description: "Reference name"
    }
    gvcf_output: {
      description: "Output GVCF?"
    }
    parabricks_version: {
      description: "Parabricks version"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File? regions_bed
    File ref_fasta
    File ref_index
    String ref_name
    Boolean gvcf_output = true
    String parabricks_version = "4.7.0-1"
    RuntimeAttributes default_runtime_attributes
  }

  String docker_image = if (default_runtime_attributes.backend == "AWS-OMICS")
    then default_runtime_attributes.container_registry
    else "nvcr.io/nvidia/clara" + "/clara-parabricks:~{parabricks_version}"

  call run_parabricks_deepvariant { input:
    sample_id = sample_id,
    aligned_bam = aligned_bam,
    aligned_bam_index = aligned_bam_index,
    regions_bed = regions_bed,
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    ref_name = ref_name,
    gvcf_output = gvcf_output,
    docker_image = docker_image,
    runtime_attributes = default_runtime_attributes
  }

  call postprocess_vcf { input:
    vcf = run_parabricks_deepvariant.vcf,
    gvcf = run_parabricks_deepvariant.gvcf,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File vcf = postprocess_vcf.vcf_gz
    File vcf_index = postprocess_vcf.vcf_index
    File? gvcf = postprocess_vcf.gvcf_gz
    File? gvcf_index = postprocess_vcf.gvcf_index
  }
}

task run_parabricks_deepvariant {
  meta {
    description: "Call variants from aligned HiFi reads using Parabricks DeepVariant"
    outputs: {
      vcf: {
        description: "Small variant VCF"
      },
      gvcf: {
        description: "Small variant GVCF"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    aligned_bam: {
      description: "Aligned BAM"
    }
    aligned_bam_index: {
      description: "Aligned BAM index"
    }
    regions_bed: {
      description: "Regions BED"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    ref_name: {
      description: "Reference name"
    }
    gvcf_output: {
      description: "Output gVCF?"
    }
    docker_image: {
      description: "Docker image URL"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File? regions_bed
    File ref_fasta
    File ref_index
    String ref_name
    Boolean gvcf_output
    String docker_image
    RuntimeAttributes runtime_attributes
  }

  Int threads = 48
  Int numGPUs = 4
  Int mem_gb = threads * 4
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    export TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456

    /usr/local/parabricks/pbrun deepvariant \
      --num-gpus ~{numGPUs} \
      --run-partition \
      --mode pacbio \
      --ref "~{ref_fasta}" \
      ~{if defined(regions_bed)
        then "--interval-file '" + regions_bed + "'"
        else ""
      } \
      ~{if (gvcf_output)
        then "--gvcf"
        else ""
      } \
      --in-bam "~{aligned_bam}" \
      --out-variants "~{sample_id}.~{ref_name}.small_variants~{if (gvcf_output)
        then ".g"
        else ""
      }.vcf"
  >>>

  output {
    File vcf = "~{sample_id}.~{ref_name}.small_variants.vcf"
    File? gvcf = "~{sample_id}.~{ref_name}.small_variants.g.vcf"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    gpuCount: numGPUs
    gpuType: runtime_attributes.gpuType
    acceleratorCount: numGPUs  # !UnknownRuntimeKey
    acceleratorType: runtime_attributes.gpuType  # !UnknownRuntimeKey
    nvidiaDriverVersion: "535.230.02"  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task postprocess_vcf {
  meta {
    description: "Filter, compress, and index VCFs"
    outputs: {
      vcf_gz: {
        description: "Compressed VCF"
      },
      vcf_index: {
        description: "VCF index"
      },
      gvcf_gz: {
        description: "Compressed GVCF"
      },
      gvcf_index: {
        description: "GVCF index"
      }
    }
  }

  parameter_meta {
    vcf: {
      description: "VCF"
    }
    gvcf: {
      description: "GVCF"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File vcf
    File? gvcf
    RuntimeAttributes runtime_attributes
  }

  Int threads = 8
  Int mem_gb = 16
  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  String out_prefix = basename(vcf, ".vcf")

  command <<<
    set -euo pipefail

    bgzip --version
    bcftools --version

    if [[ "~{defined(gvcf)}" == "true" ]]; then
      bgzip --stdout \
        ~{if threads > 1
          then "--threads '" + (threads - 1) + "'"
          else ""
        } \
        "~{gvcf}" \
      > "~{out_prefix}.g.vcf.gz"
      bcftools index --tbi --force \
        ~{if threads > 1
          then "--threads '" + (threads - 1) + "'"
          else ""
        } \
        "~{out_prefix}.g.vcf.gz"
    fi

    # Filter VCF to remove uncalled sites
    bcftools view \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --exclude-uncalled \
      --output-type z \
      --output-file "~{out_prefix}.vcf.gz" \
      "~{vcf}"
    bcftools index --tbi --force \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      "~{out_prefix}.vcf.gz"
  >>>

  output {
    File vcf_gz = "~{out_prefix}.vcf.gz"
    File vcf_index = "~{out_prefix}.vcf.gz.tbi"
    File? gvcf_gz = "~{out_prefix}.g.vcf.gz"
    File? gvcf_index = "~{out_prefix}.g.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"  # pb_wdl_base:build2
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}
