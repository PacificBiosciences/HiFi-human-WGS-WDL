# Reference Map File Specification

| Type | Key | Description | Notes |
| ---- | --- | ----------- | ----- |
| String | name | Short name for reference | Alphanumeric characters, underscores, and dashes only.  Will be used in file names. |
| File | fasta | Reference genome FASTA |  |
| File | fasta_index | Reference genome FASTA index |  |
| File | pbsv_splits | Regions for pbsv parallelization | [below](#pbsv_splits) |
| File | pbsv_tandem_repeat_bed | Tandem Repeat BED used by PBSV to normalize SVs within TRs | [link](https://github.com/PacificBiosciences/pbsv/tree/master/annotations) |
| File | trgt_tandem_repeat_bed | Tandem Repeat catalog (BED) for TRGT genotyping | [link](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md) |
| File | sawfish_exclude_bed | Regions to be excluded for Sawfish CNV calls in gzipped BED format | [link](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md#cnv-excluded-regions) |
| File | sawfish_exclude_bed_index | BED index | [link](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md#cnv-excluded-regions) |
| File | sawfish_expected_bed_male | Expected allosome copy number BED for XY samples | [link](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md#expected-copy-number) |
| File | sawfish_expected_bed_female | Expected allosome copy number BED for XX samples | [link](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md#expected-copy-number) |
| File | pharmcat_positions_vcf | PharmCAT positions VCF |  |
| File | pharmcat_positions_vcf_index | PharmCAT positions VCF index |  |
