version 1.0

import "../structs.wdl"

task write_ped_phrank {
  meta {
    description: "Create a tab-delimited PLINK pedigree (PED) format and calculate Phrank phenotype match score for every gene."
  }

  parameter_meta {
    id: {
      name: "Family ID if passing a Family JSON struct. Sample ID if passing a single sample."
    }
    family_json: {
      name: "Family JSON struct"
    }
    phenotypes: {
      name: "Comma delimited string of HPO terms for phenotypes."
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pedigree: {
      name: "PLINK pedigree (PED) format"
    }
    phrank_lookup: {
      name: "Phrank scores: <gene_symbol><TAB><phrank_score>"
    }
  }

  input {
    String id

    Array[Map[String, String]] family_json

    String phenotypes

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 1
  Int mem_gb    = 2
  Int disk_size = 1

  command <<<
    set -euo pipefail

    json2ped.py ~{write_json(family_json)} ~{id} > ~{id}.ped

    cat ~{id}.ped

    # clean up common potential issues with phenotypes string
    PHENOTYPES="~{phenotypes}"
    PHENOTYPES="${PHENOTYPES// /}"   # Remove spaces
    PHENOTYPES="${PHENOTYPES//;/,}"  # Replace semicolons with commas

    # ENV HPO_TERMS_TSV "/opt/data/hpo/hpoTerms.txt"
    # ENV HPO_DAG_TSV "/opt/data/hpo/hpoDag.txt"
    # ENV ENSEMBL_TO_HPO_TSV "/opt/data/hpo/ensembl.hpoPhenotype.tsv"
    # ENV ENSEMBL_TO_HGNC "/opt/data/genes/ensembl.hgncSymbol.tsv"

    calculate_phrank.py \
      "${HPO_TERMS_TSV}" \
      "${HPO_DAG_TSV}" \
      "${ENSEMBL_TO_HPO_TSV}" \
      "${ENSEMBL_TO_HGNC}" \
      "${PHENOTYPES}" \
      ~{id}_phrank.tsv
  >>>

  output {
    File pedigree      = "~{id}.ped"
    File phrank_lookup = "~{id}_phrank.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/wgs_tertiary@sha256:6925e50c6850d8886a9cbe0d984ec8eaea205cca280309820688b60dfc65b1d4"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}