# Reference Map File Specification

| Type | Key | Description | Notes |
| ---- | --- | ----------- | ----- |
| String | name | Short name for reference | Alphanumeric characters, underscores, and dashes only.  Will be used in file names. |
| File | fasta | Reference genome FASTA |  |
| File | fasta_index | Reference genome FASTA index |  |
| File | pbsv_splits | Regions for pbsv parallelization | [below](#pbsv_splits) |
| File | pbsv_tandem_repeat_bed | Tandem Repeat BED used by PBSV to normalize SVs within TRs | [link](https://github.com/PacificBiosciences/pbsv/tree/master/annotations) |
| File | trgt_tandem_repeat_bed | Tandem Repeat catalog (BED) for TRGT genotyping | [link](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md) |
| File | hificnv_exclude_bed | Regions to be excluded by HIFICNV in gzipped BED format | [link](https://github.com/PacificBiosciences/HiFiCNV/blob/main/docs/aux_data.md) |
| File | hificnv_exclude_bed_index | BED index | [link](https://github.com/PacificBiosciences/HiFiCNV/blob/main/docs/aux_data.md) |
| File | hificnv_expected_bed_male | Expected allosome copy number BED for XY samples | [link](https://github.com/PacificBiosciences/HiFiCNV/blob/main/docs/aux_data.md) |
| File | hificnv_expected_bed_female | Expected allosome copy number BED for XX samples | [link](https://github.com/PacificBiosciences/HiFiCNV/blob/main/docs/aux_data.md) |
| File | pharmcat_positions_vcf | PharmCAT positions VCF |  |
| File | pharmcat_positions_vcf_index | PharmCAT positions VCF index |  |

## pbsv_splits

The `pbsv_splits` file is a JSON array of arrays of strings. Each inner array contains one or more chromosome names such that each inner array is of roughly equal size in base pairs. The inner arrays are processed in parallel.  For example:

```json
[
  ...
    [
        "chr10",
        "chr11"
    ],
    [
        "chr12",
        "chr13"
    ],
    [
        "chr14",
        "chr15"
    ],
  ...
]
```
