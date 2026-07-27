version 1.0

import "../structs.wdl"

task glnexus {
  meta {
    description: "Merge gVCFs from multiple samples and joint call genotypes with GLnexus"
    outputs: {
      vcf: {
        description: "Joint-called multi-sample VCF"
      },
      vcf_index: {
        description: "Index for joint-called multi-sample VCF"
      }
    }
  }

  parameter_meta {
    cohort_id: {
      description: "Cohort ID"
    }
    gvcfs: {
      description: "GVCFs"
    }
    gvcf_indices: {
      description: "GVCF indices"
    }
    ref_name: {
      description: "Reference name"
    }
    regions_bed: {
      description: "Regions BED"
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
    String cohort_id
    Array[File] gvcfs
    Array[File] gvcf_indices
    String ref_name
    File? regions_bed
    Int threads = 32
    Int mem_gb = 60
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(gvcfs, "GB") * 2 + 100)

  command <<<
    set -euo pipefail

    # we use a custom config file based on DeepVariant_unfiltered, but with required_dp: 1
    cat << EOF > config.yml
    unifier_config:
      min_AQ1: 0
      min_AQ2: 0
      min_GQ: 0
      monoallelic_sites_for_lost_alleles: true
      max_alleles_per_site: 32
    genotyper_config:
      required_dp: 1
      revise_genotypes: false
      allow_partial_data: true
      more_PL: true
      trim_uncalled_alleles: true
      liftover_fields:
        - orig_names: [MIN_DP, DP]
          name: DP
          description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
          type: int
          combi_method: min
          number: basic
          count: 1
          ignore_non_variants: true
        - orig_names: [AD]
          name: AD
          description: '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
          type: int
          number: alleles
          combi_method: min
          default_type: zero
          count: 0
        - orig_names: [GQ]
          name: GQ
          description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
          type: int
          number: basic
          combi_method: min
          count: 1
          ignore_non_variants: true
        - orig_names: [PL]
          name: PL
          description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
          type: int
          number: genotype
          combi_method: missing
          count: 0
          ignore_non_variants: true
    EOF

    GVCFS=()
    for i in ~{sep=" " gvcfs}; do
      ln --symbolic --verbose "${i}" .
      # shellcheck disable=SC2086
      GVCFS+=("$(basename ${i})")
    done
    for i in ~{sep=" " gvcf_indices}; do
      ln --symbolic --verbose "${i}" .
    done

    ~{if defined(regions_bed)
      then "ln --symbolic --verbose '" + regions_bed + "' ."
      else ""
    }

    # glneux_cli has no version option
    glnexus_cli --help 2>&1 | grep -Eo 'glnexus_cli release v[0-9a-f.-]+'
    bcftools --version

    # shellcheck disable=SC2068
    glnexus_cli \
      --threads ~{threads} \
      --mem-gbytes ~{mem_gb} \
      --dir "~{cohort_id}.~{ref_name}.GLnexus.DB" \
      --config ./config.yml \
      ~{if defined(regions_bed)
        then "--bed '" + basename(select_first([
          regions_bed
        ])) + "'"
        else ""
      } \
      ${GVCFS[@]} \
    > "~{cohort_id}.~{ref_name}.small_variants.bcf"
    bcftools view \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --output-type z \
      --output-file "~{cohort_id}.~{ref_name}.small_variants.vcf.gz" \
      "~{cohort_id}.~{ref_name}.small_variants.bcf"
    bcftools index \
      ~{if threads > 1
        then "--threads '" + (threads - 1) + "'"
        else ""
      } \
      --tbi \
      "~{cohort_id}.~{ref_name}.small_variants.vcf.gz"
    # cleanup intermediate files
    rm --recursive --force --verbose "~{cohort_id}.~{ref_name}.GLnexus.DB" "~{cohort_id}.~{ref_name}.small_variants.bcf"
  >>>

  output {
    File vcf = "~{cohort_id}.~{ref_name}.small_variants.vcf.gz"
    File vcf_index = "~{cohort_id}.~{ref_name}.small_variants.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/glnexus@sha256:ce6fecf59dddc6089a8100b31c29c1e6ed50a0cf123da9f2bc589ee4b0c69c8e"  # 1.4.3
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

