version 1.0

import "../../structs.wdl"

workflow deepvariant {
  meta {
    description: "Call variants from aligned HiFi reads using DeepVariant"
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
    aligned_bams: {
      description: "Aligned BAMs"
    }
    aligned_bam_indices: {
      description: "Aligned BAM indices"
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
    output_deepvariant_phasing: {
      description: "Output phased variants using DeepVariant continuous phasing information?"
    }
    gvcf_output: {
      description: "Output GVCF?"
    }
    deepvariant_version: {
      description: "DeepVariant version"
    }
    gpu: {
      description: "Use GPU for DeepVariant"
    }
    max_nproc: {
      description: "Max CPUs available on a single node in the target cluster -- drives call_variants_cpu's threads (min(max_nproc, 64)); does not affect call_variants_gpu or postprocess_variants, which don't scale with node size"
    }
    default_runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    File? regions_bed
    File ref_fasta
    File ref_index
    String ref_name
    Boolean output_deepvariant_phasing = false
    Boolean gvcf_output = true
    String deepvariant_version = "1.10.0"
    Boolean gpu
    Int max_nproc = 64
    RuntimeAttributes default_runtime_attributes
  }

  Int total_deepvariant_tasks = 64
  Int num_shards = 8
  Int tasks_per_shard = total_deepvariant_tasks / num_shards

  String docker_image = "google/deepvariant:~{deepvariant_version}"

  scatter (shard_index in range(num_shards)) {
    Int task_start_index = shard_index * tasks_per_shard

    call deepvariant_make_examples { input:
      sample_id = sample_id,
      aligned_bams = aligned_bams,
      aligned_bam_indices = aligned_bam_indices,
      regions_bed = regions_bed,
      ref_fasta = ref_fasta,
      ref_index = ref_index,
      task_start_index = task_start_index,
      tasks_per_shard = tasks_per_shard,
      total_deepvariant_tasks = total_deepvariant_tasks,
      output_deepvariant_phasing = output_deepvariant_phasing,
      gvcf_output = gvcf_output,
      docker_image = docker_image,
      runtime_attributes = default_runtime_attributes
    }
  }

  if (!gpu) {
    call deepvariant_call_variants_cpu { input:
      sample_id = sample_id,
      ref_name = ref_name,
      example_tfrecord_tars = deepvariant_make_examples.example_tfrecord_tar,
      total_deepvariant_tasks = total_deepvariant_tasks,
      docker_image = docker_image,
      threads = if (max_nproc < 64)
        then max_nproc
        else 64
      ,
      runtime_attributes = default_runtime_attributes
    }
  }

  if (gpu) {
    call deepvariant_call_variants_gpu { input:
      sample_id = sample_id,
      ref_name = ref_name,
      example_tfrecord_tars = deepvariant_make_examples.example_tfrecord_tar,
      total_deepvariant_tasks = total_deepvariant_tasks,
      docker_image = docker_image + "-gpu",
      runtime_attributes = default_runtime_attributes
    }
  }

  call deepvariant_postprocess_variants { input:
    sample_id = sample_id,
    tfrecords_tar = select_first([
      deepvariant_call_variants_gpu.tfrecords_tar,
      deepvariant_call_variants_cpu.tfrecords_tar
    ]),
    example_tfrecord_tars = deepvariant_make_examples.example_tfrecord_tar,
    nonvariant_site_tfrecord_tars = select_all(deepvariant_make_examples.nonvariant_site_tfrecord_tar),
    read_phasing_tars = select_all(deepvariant_make_examples.read_phasing_tar),
    ref_fasta = ref_fasta,
    ref_index = ref_index,
    ref_name = ref_name,
    total_deepvariant_tasks = total_deepvariant_tasks,
    docker_image = docker_image,
    runtime_attributes = default_runtime_attributes
  }

  output {
    File vcf = deepvariant_postprocess_variants.vcf
    File vcf_index = deepvariant_postprocess_variants.vcf_index
    File? gvcf = deepvariant_postprocess_variants.gvcf
    File? gvcf_index = deepvariant_postprocess_variants.gvcf_index
  }
}

task deepvariant_make_examples {
  meta {
    description: "Run DeepVariant make_examples step"
    outputs: {
      example_tfrecord_tar: {
        description: "Example TFRecord tar"
      },
      nonvariant_site_tfrecord_tar: {
        description: "Nonvariant site TFRecord tar"
      },
      read_phasing_tar: {
        description: "Read phasing tar"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    aligned_bams: {
      description: "Aligned BAM"
    }
    aligned_bam_indices: {
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
    task_start_index: {
      description: "Task start index"
    }
    tasks_per_shard: {
      description: "Tasks per shard"
    }
    total_deepvariant_tasks: {
      description: "Total DeepVariant tasks"
    }
    output_deepvariant_phasing: {
      description: "Phase variants with DeepVariant continuous phasing information"
    }
    gvcf_output: {
      description: "Output GVCF nonvariant site TFRecords"
    }
    docker_image: {
      description: "Docker image URL"
    }
    threads_per_task: {
      description: "Number of threads to use per task"
    }
    mem_gb_per_task: {
      description: "Memory to use in GiB per task"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    File? regions_bed
    File ref_fasta
    File ref_index
    Int task_start_index
    Int tasks_per_shard
    Int total_deepvariant_tasks
    Boolean output_deepvariant_phasing
    Boolean gvcf_output
    String docker_image
    Int threads_per_task = 1
    Int mem_gb_per_task = 4
    RuntimeAttributes runtime_attributes
  }

  Int task_end_index = task_start_index + tasks_per_shard - 1

  Int disk_size = ceil(size(aligned_bams, "GB") * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    BAMS=()
    for i in ~{sep=" " aligned_bams}; do
      ln --symbolic --verbose "${i}" .
      # shellcheck disable=SC2086
      BAMS+=("$(basename ${i})")
    done
    for i in ~{sep=" " aligned_bam_indices}; do
      ln --symbolic --verbose "${i}" .
    done

    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .
    ~{if defined(regions_bed)
      then "ln --symbolic --verbose '" + regions_bed + "' ."
      else ""
    }

    mkdir example_tfrecords nonvariant_site_tfrecords read_phasing

    # shellcheck disable=SC2068,SC2046,SC2086
    seq ~{task_start_index} ~{task_end_index} \
    | parallel \
      --jobs ~{tasks_per_shard} \
      --halt 2 \
      /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref "~{basename(ref_fasta)}" \
        ~{if defined(regions_bed)
          then "--regions ~{basename(select_first([
            regions_bed
          ]))}"
          else ""
        } \
        --reads $(IFS=,; echo "${BAMS[*]}") \
        --examples "example_tfrecords/make_examples.tfrecord@~{total_deepvariant_tasks}.gz" \
        ~{if (gvcf_output)
          then "--gvcf 'nonvariant_site_tfrecords/gvcf.tfrecord@~{total_deepvariant_tasks}.gz'"
          else ""
        } \
        ~{if (output_deepvariant_phasing)
          then "--output_phase_info=true"
          else ""
        } \
        ~{if (output_deepvariant_phasing)
          then "--output_local_read_phasing='read_phasing/read_phasing_debug@~{total_deepvariant_tasks}.tsv'"
          else ""
        } \
        --checkpoint /opt/models/pacbio \
        --task {}

    tar --gzip --create --verbose --file "~{sample_id}.~{task_start_index}.example_tfrecords.tar.gz" example_tfrecords \
    && rm --recursive --force --verbose example_tfrecords
    ~{if (gvcf_output)
      then "tar --gzip --create --verbose --file ~{sample_id}.~{task_start_index}.nonvariant_site_tfrecords.tar.gz nonvariant_site_tfrecords && rm --recursive --force --verbose nonvariant_site_tfrecords"
      else ""
    }
    ~{if (output_deepvariant_phasing)
      then "tar --gzip --create --verbose --file ~{sample_id}.~{task_start_index}.read_phasing.tar.gz read_phasing && rm --recursive --force --verbose read_phasing"
      else ""
    }
  >>>

  output {
    File example_tfrecord_tar = "~{sample_id}.~{task_start_index}.example_tfrecords.tar.gz"
    File? nonvariant_site_tfrecord_tar = "~{sample_id}.~{task_start_index}.nonvariant_site_tfrecords.tar.gz"
    File? read_phasing_tar = "~{sample_id}.~{task_start_index}.read_phasing.tar.gz"
  }

  runtime {
    docker: docker_image
    cpu: tasks_per_shard * threads_per_task
    memory: "~{tasks_per_shard * mem_gb_per_task} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task deepvariant_call_variants_cpu {
  meta {
    description: "Run DeepVariant call_variants step"
    outputs: {
      tfrecords_tar: {
        description: "TFRecords tar"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    ref_name: {
      description: "Reference name"
    }
    example_tfrecord_tars: {
      description: "Example TFRecord tars"
    }
    total_deepvariant_tasks: {
      description: "Total DeepVariant tasks"
    }
    docker_image: {
      description: "Docker image URL"
    }
    threads: {
      description: "Number of threads to use"
    }
    writer_threads: {
      description: "Number of writer threads to use"
    }
    batch_size: {
      description: "call_variants --batch_size, DV default is 1024."
    }
    mem_gb: {
      description: "Memory to use in GiB. Default is a formula of batch_size, override to set explicitly"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String ref_name
    Array[File] example_tfrecord_tars
    Int total_deepvariant_tasks
    String docker_image
    Int threads = 64
    Int writer_threads = 8
    Int batch_size = 1024
    # mem_gb is batch_size-driven, not threads-driven: true-peak-fit + 1.5x
    # margin, rounded to the nearest 4 GiB. 28 GiB at the default batch_size.
    Int mem_gb = 4 * ceil(1.5 * (12.50 + 0.00502 * batch_size) / 4)
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(example_tfrecord_tars, "GB") * 2 + 100)

  command <<<
    set -euo pipefail

    while read -r tfrecord_tar || [[ -n "${tfrecord_tar}" ]]; do
      tar --no-same-owner --gzip --extract --verbose --file "${tfrecord_tar}"
    done < "~{write_lines(example_tfrecord_tars)}"

    /opt/deepvariant/bin/call_variants \
      --writer_threads ~{writer_threads} \
      --batch_size ~{batch_size} \
      --outfile call_variants_output.tfrecord.gz \
      --examples "example_tfrecords/make_examples.tfrecord@~{total_deepvariant_tasks}.gz" \
      --checkpoint "/opt/models/pacbio"

    tar --gzip --create --verbose --file "~{sample_id}.~{ref_name}.call_variants_output.tar.gz" call_variants_output*.tfrecord.gz \
    && rm --verbose call_variants_output*.tfrecord.gz \
    && rm --recursive --force --verbose example_tfrecords
  >>>

  output {
    File tfrecords_tar = "~{sample_id}.~{ref_name}.call_variants_output.tar.gz"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task deepvariant_call_variants_gpu {
  meta {
    description: "Run DeepVariant call_variants step"
    outputs: {
      tfrecords_tar: {
        description: "TFRecords tar"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    ref_name: {
      description: "Reference name"
    }
    example_tfrecord_tars: {
      description: "Example TFRecord tars"
    }
    total_deepvariant_tasks: {
      description: "Total DeepVariant tasks"
    }
    docker_image: {
      description: "Docker image URL"
    }
    threads: {
      description: "Number of threads to use"
    }
    writer_threads: {
      description: "Number of writer threads to use"
    }
    mem_gb: {
      description: "Memory to use in GiB."
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String ref_name
    Array[File] example_tfrecord_tars
    Int total_deepvariant_tasks
    String docker_image
    Int threads = 8
    Int writer_threads = 4
    Int mem_gb = 44
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(example_tfrecord_tars, "GB") * 2 + 100)
  Int gpuCount = 1

  command <<<
    set -euo pipefail

    while read -r tfrecord_tar || [[ -n "${tfrecord_tar}" ]]; do
      tar --no-same-owner --gzip --extract --verbose --file "${tfrecord_tar}"
    done < "~{write_lines(example_tfrecord_tars)}"

    /opt/deepvariant/bin/call_variants \
      --writer_threads ~{writer_threads} \
      --outfile call_variants_output.tfrecord.gz \
      --examples "example_tfrecords/make_examples.tfrecord@~{total_deepvariant_tasks}.gz" \
      --checkpoint "/opt/models/pacbio"

    tar --gzip --create --verbose --file "~{sample_id}.~{ref_name}.call_variants_output.tar.gz" call_variants_output*.tfrecord.gz \
    && rm --verbose call_variants_output*.tfrecord.gz \
    && rm --recursive --force --verbose example_tfrecords
  >>>

  output {
    File tfrecords_tar = "~{sample_id}.~{ref_name}.call_variants_output.tar.gz"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 30  # !UnknownRuntimeKey
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    gpuCount: gpuCount
    gpuType: runtime_attributes.gpuType
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

task deepvariant_postprocess_variants {
  meta {
    description: "Run DeepVariant postprocess_variants step"
    outputs: {
      vcf: {
        description: "VCF"
      },
      vcf_index: {
        description: "VCF index"
      },
      gvcf: {
        description: "GVCF"
      },
      gvcf_index: {
        description: "GVCF index"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    tfrecords_tar: {
      description: "TFRecords tar"
    }
    example_tfrecord_tars: {
      description: "Example TFRecord tars"
    }
    nonvariant_site_tfrecord_tars: {
      description: "Nonvariant site TFRecord tars"
    }
    read_phasing_tars: {
      description: "Read phasing tar"
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
    total_deepvariant_tasks: {
      description: "Total DeepVariant tasks"
    }
    docker_image: {
      description: "Docker image URL"
    }
    threads: {
      description: "Number of threads to use"
    }
    mem_gb: {
      description: "Memory to use in GiB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    File tfrecords_tar
    Array[File] example_tfrecord_tars
    Array[File?] nonvariant_site_tfrecord_tars
    Array[File?] read_phasing_tars
    File ref_fasta
    File ref_index
    String ref_name
    Int total_deepvariant_tasks
    String docker_image
    Int threads = 8
    Int mem_gb = 48
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil((size(ref_fasta, "GB") + size(tfrecords_tar, "GB") + size(example_tfrecord_tars, "GB") + (if length(nonvariant_site_tfrecord_tars) > 0
    then size(select_all(nonvariant_site_tfrecord_tars), "GB")
    else 0
  ) + (if length(read_phasing_tars) > 0
    then size(select_all(read_phasing_tars), "GB")
    else 0
  )) * 2 + 20)

  command <<<
    set -euo pipefail

    tar --no-same-owner --gzip --extract --verbose --file "~{tfrecords_tar}"

    while read -r tfrecord_tar || [[ -n "${tfrecord_tar}" ]]; do
      tar --no-same-owner --gzip --extract --verbose --file "${tfrecord_tar}"
    done < "~{write_lines(example_tfrecord_tars)}"

    if [ "~{length(read_phasing_tars) > 0}" == "true" ]; then
      while read -r read_phasing_tar || [[ -n "${read_phasing_tar}" ]]; do
        tar --no-same-owner --gzip --extract --verbose --file "${read_phasing_tar}"
      done < "~{write_lines(select_all(read_phasing_tars))}"
    fi

    if [ "~{length(nonvariant_site_tfrecord_tars) > 0}" == "true" ]; then
      while read -r nonvariant_site_tfrecord_tar || [[ -n "${nonvariant_site_tfrecord_tar}" ]]; do
        tar --no-same-owner --gzip --extract --verbose --file "${nonvariant_site_tfrecord_tar}"
      done < "~{write_lines(select_all(nonvariant_site_tfrecord_tars))}"
    fi

    ln --symbolic --verbose "~{ref_fasta}" .
    ln --symbolic --verbose "~{ref_index}" .

    # shellcheck disable=SC2086
    /opt/deepvariant/bin/postprocess_variants \
      --cpus ~{threads} \
      --vcf_stats_report=false \
      --ref "~{basename(ref_fasta)}" \
      --infile call_variants_output.tfrecord.gz \
      --small_model_cvo_records "example_tfrecords/make_examples_call_variant_outputs.tfrecord@~{total_deepvariant_tasks}.gz" \
      ~{if (length(nonvariant_site_tfrecord_tars) > 0)
        then "--nonvariant_site_tfrecord_path nonvariant_site_tfrecords/gvcf.tfrecord@~{total_deepvariant_tasks}.gz --gvcf_outfile ~{sample_id}.~{ref_name}.small_variants.g.vcf.gz"
        else ""
      } \
      ~{if (length(read_phasing_tars) > 0)
        then "--phased_reads_input_path=read_phasing/read_phasing_debug@~{total_deepvariant_tasks}.tsv"
        else ""
      } \
      --outfile "~{sample_id}.~{ref_name}.small_variants.vcf.gz"

    # Filter for only PASS variants
    bcftools view \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --exclude-uncalled \
      --output-type z \
      --output-file "~{sample_id}.~{ref_name}.small_variants.passing.vcf.gz" \
      "~{sample_id}.~{ref_name}.small_variants.vcf.gz"

    mv --verbose "~{sample_id}.~{ref_name}.small_variants.passing.vcf.gz" "~{sample_id}.~{ref_name}.small_variants.vcf.gz"
    bcftools index --tbi --force \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      "~{sample_id}.~{ref_name}.small_variants.vcf.gz"

    rm --verbose call_variants_output*.tfrecord.gz \
    && rm --recursive --force --verbose example_tfrecords nonvariant_site_tfrecords read_phasing
  >>>

  output {
    File vcf = "~{sample_id}.~{ref_name}.small_variants.vcf.gz"
    File vcf_index = "~{sample_id}.~{ref_name}.small_variants.vcf.gz.tbi"
    File? gvcf = "~{sample_id}.~{ref_name}.small_variants.g.vcf.gz"
    File? gvcf_index = "~{sample_id}.~{ref_name}.small_variants.g.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{mem_gb} GiB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}

