version 1.0

import "../structs.wdl"
import "../../../humanwgs_structs.wdl"

task write_ped_phrank {
  meta {
    description: "Create a tab-delimited PLINK pedigree (PED) format and calculate Phrank phenotype match score for every gene."
  }

  parameter_meta {
    id: {
      name: "Family ID if passing a Family struct. Sample ID if passing a single sample."
    }
    family: {
      name: "Family struct (required if passing a Family)."
    }
    sex: {
      name: "Sample sex (required if passing a single sample)."
    }
    phenotypes: {
      name: "Comma delimited string of HPO terms for phenotypes."
    }
    disk_size: {
      name: "Disk size in GB"
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

    Family? family
    String? sex

    String phenotypes

    Int disk_size = 1

    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb  = 2

  command <<<
    set -euo pipefail

    if ~{defined(family)}; then
      echo "Family struct provided. Converting to PED format."
    json2ped.py <(jq '.[0]' ~{write_json(select_all([family]))}) > ~{id}.ped
    else
      echo "Family struct not provided. Creating single sample PED file."
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
    docker: "~{runtime_attributes.container_registry}/wgs_tertiary@sha256:410597030e0c85cf16eb27a877d260e7e2824747f5e8b05566a1aaa729d71136"
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