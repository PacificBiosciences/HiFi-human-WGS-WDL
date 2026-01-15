version 1.0

import "../structs.wdl"

task write_phrank {
  meta {
    description: "Calculate Phrank phenotype match score for every gene."
  }

  parameter_meta {
    phenotypes: {
      name: "Comma delimited string of HPO terms for phenotypes."
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    phrank_lookup: {
      name: "Phrank scores: <gene_symbol><TAB><phrank_score>"
    }
  }

  input {
    String phenotypes

    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb  = 2
  Int disk_size = 1

  command <<<
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
      phrank.tsv
  >>>

  output {
    File phrank_lookup = "phrank.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/wgs_tertiary@sha256:410597030e0c85cf16eb27a877d260e7e2824747f5e8b05566a1aaa729d71136"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}