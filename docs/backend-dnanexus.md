# PacBio Human WGS Variant Pipeline on DNAnexus

The PacBio Human WGS Variant Pipeline is an analysis workflow for PacBio HiFi human whole genome sequencing data, with joint calling for related samples.

Templates and instructions for how to submit the `family` input on the DNAnexus platform are provided in the [Example JSON Documents](#example-json-documents) and [Submitting to DNAnexus](#submitting-to-dnanexus) sections below.  

## Inputs

The workflow has the following inputs:

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Struct | family | Family struct describing samples, relationships, and unaligned BAM paths. | See below for more information |
| String | phenotypes | [Human Phenotype Ontology](https://hpo.jax.org/) (HPO) phenotypes associated with the affected proband. | For example, if the proband has seizures and hypotonia, then the `phenotypes` string might be `"HP:0001250,HP:0001252"`. |
| File | trgt_tandem_repeat_bed | BED file containing repeat coordinates and information about the repeat structure | The default file should be sufficient for most use cases |
| Integer | glnexus_mem_gb | Override GLnexus memory request (GB) | Optional. Should only be specified if GLnexus step fails. |
| Integer | pbsv_call_mem_gb | Override PBSV call memory request (GB) | Optional. Should only be specified if PBSV step fails. |
| Boolean | run_tertiary | Whether to run tertiary analysis for small variants and structural variants | Default: `true` |

## Family Struct Syntax

The `Family` input for the HiFi-human-WGS-WDL workflow is a JSON document that contains the samples for the family. The same struct is used for a single sample or trio, with the single sample case having only one entry in the `samples` array.

### Structs and Field Descriptions

#### Family Struct

The `Family` struct contains the samples for the family. The struct has the following fields:

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | family_id | Unique identifier for the family | Alphanumeric characters, periods, dashes, and underscores are allowed. |
| Array\[[Sample](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/workflows/humanwgs_structs.wdl#L3)\] | samples | Sample struct with sample specific data and metadata. | [below](#sample-struct) |

#### Sample Struct

A `Sample` struct contains sample specific data and metadata. The struct has the following fields:

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | sample_id | Unique identifier for the sample | Alphanumeric characters, periods, dashes, and underscores are allowed. |
| String? | sex | Sample sex<br/>`["MALE", "FEMALE"]` | Optional field used by HiFiCNV and TRGT for genotyping. Allosome karyotype will default to XX unless sex is specified as `"MALE"`.  Used for tertiary analysis X-linked inheritance filtering. |
| Boolean | affected | Affected status | If set to `true`, sample is described as being affected by all HPO terms in `phenotypes`.<br/>If set to `false`, sample is described as not being affected by all HPO terms in `phenotypes`. |
| Array\[File\] | hifi_reads | Array of [DNAnexus links](https://documentation.dnanexus.com/user/projects/path-resolution#dnanexus-links) for HiFi reads in unaligned BAM format. | |
| String? | father_id | sample_id of father (optional) | |
| String? | mother_id | sample_id of mother (optional) | |

### Example JSON Documents

#### Specifying HiFi BAM inputs

When specifying DNAnexus `hifi_reads` files, format them as `{"$dnanexus_link": {"id": "file-xxxx", "project": "project-xxxx"}}` as is in the examples below. One or more `hifi_reads` files can be specified for each sample.  The single sample example has multiple files, and the trio example has one file per sample.

#### Example Single Sample JSON

In this example, the optional `sex` field is not specified, so tools will default to XX for the allosome karyotype.

```json
{
  "family_id": "EXAMPLE-singleton",
  "samples": [
    {
      "sample_id": "EXAMPLE",
      "hifi_reads": [
        {
          "$dnanexus_link": {
            "id": "file-xxxx",
            "project": "project-xxxx"
          }
        },
        {
          "$dnanexus_link": {
            "id": "file-xxxx",
            "project": "project-xxxx"
          }
        }
      ],
      "affected": true
    }
  ]
}
```

#### Example Trio JSON

```json
{
  "family": {
    "family_id": "AJTRIO",
    "samples": [
      {
        "sample_id": "HG002",
        "hifi_reads": [
          {
            "$dnanexus_link": {
              "id": "file-xxxx",
              "project": "project-xxxx"
            }
          }
        ],
        "affected": true,
        "sex": "MALE",
        "father_id": "HG003",
        "mother_id": "HG004"
      },
      {
        "sample_id": "HG003",
        "hifi_reads": [
          {
            "$dnanexus_link": {
              "id": "file-xxxx",
              "project": "project-xxxx"
            }
          }
        ],
        "affected": false,
        "sex": "MALE"
      },
      {
        "sample_id": "HG004",
        "hifi_reads": [
            {
            "$dnanexus_link": {
              "id": "file-xxxx",
              "project": "project-xxxx"
            }
          }
        ],
        "affected": false,
        "sex": "FEMALE"
      }
    ]
  }
}
```

### Submitting to DNAnexus

When submitting files to DNAnexus, start by manually writing the JSON document for your job. The examples above can serve as a starting point. Validating the JSON for correctness can be done with an online validator such as [JSONLint](https://jsonlint.com/) or [JSONChecker](https://jsonchecker.com/).

If the job is being submitted via the DNAnexus CLI using an inputs JSON file ([documentation here](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-apps-and-applets#from-the-cli)), the JSON document should be included as the value for the `family` input parameter.

If the job is being submitted via the DNAnexus web interface, there will be two input fields for the `family` input parameter: a file array entry and a text entry. Paste the JSON document into the `family` text input field.  It is safe to paste directly in as a multi-line file as in the examples above.  In the file array entry, select all the read files that are included in the `hifi_reads` section of all included samples.  For the singleton example above, the two read files for the sample would be selected.  For the trio example above, the three read files for the three samples would be selected.

If read files included in the `family` JSON document are not selected, the workflow will fail with the error message `keys (TreeSet(affected, hifi_reads, sample_id)) have members that do not appear in struct Family`.
