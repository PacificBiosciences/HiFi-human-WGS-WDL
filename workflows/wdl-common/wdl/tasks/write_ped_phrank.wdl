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
    sex: {
      name: "Sample sex if passing a single sample."
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

    File? family_json
    String? sex

    String phenotypes

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = 20

  command <<<
    set -euo pipefail

    if ~{defined(family_json)}; then
      echo "Family JSON file provided. Converting to PED format."
      json2ped.py ~{family_json} > ~{id}.ped
    else
      echo "Family JSON file not provided. Creating single sample PED file."
      # shellcheck disable=SC2194
      case ~{if defined(sex) then sex else "."} in
        MALE | M | male | m | Male)
          SEX="1"
          ;;
        FEMALE | F | female | f | Female)
          SEX="2"
          ;;
        *)
          SEX="."
          ;;
      esac
      echo -e "~{id}\t~{id}\t.\t.\t$SEX\t2" > ~{id}.ped
    fi

    # ENV HPO_TERMS_TSV "/opt/data/hpo/hpoTerms.txt"
    # ENV HPO_DAG_TSV "/opt/data/hpo/hpoDag.txt"
    # ENV ENSEMBL_TO_HPO_TSV "/opt/data/hpo/ensembl.hpoPhenotype.tsv"
    # ENV ENSEMBL_TO_HGNC "/opt/data/genes/ensembl.hgncSymbol.tsv"

    cat ~{id}.ped

    calculate_phrank.py \
      "${HPO_TERMS_TSV}" \
      "${HPO_DAG_TSV}" \
      "${ENSEMBL_TO_HPO_TSV}" \
      "${ENSEMBL_TO_HGNC}" \
      "~{phenotypes}" \
      ~{id}_phrank.tsv
  >>>

  output {
    File pedigree      = "~{id}.ped"
    File phrank_lookup = "~{id}_phrank.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/wgs_tertiary@sha256:8fc134fdf0665e14a67bf7a8b4b63f5ae891a370a1d50c9eec2059702440a3e2"
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