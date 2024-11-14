# Tertiary Map File Specification

| Type | Key | Description | Notes |
| ---- | --- | ----------- | ----- |
| File | slivar_js | slivar functions | [link](https://raw.githubusercontent.com/brentp/slivar/91a40d582805d6607fa8a76a8fce15fd2e4be3b8/js/slivar-functions.js) |
| File | ensembl_gff | [Ensembl](https://useast.ensembl.org/index.html) GFF3 reference annotation |  |
| File | lof_lookup | Path to table of loss-of-function scores per gene |  |
| File | clinvar_lookup | Path to table of ClinVar annotations per gene |  |
| File | slivar_gnotate_files | Comma-delimited array of population dataset allele frequencies in [`slivar gnotate`](https://github.com/brentp/slivar/wiki/gnotate) format |  |
| String | slivar_gnotate_prefixes | Comma-delimieted array of prefixes to `_af`, `_nhomalt`, and `_ac` in `slivar_gnotate_files` |  |
| String (Float) [^1] | slivar_max_af | Maximum allele frequency within population for small variants |  |
| String (Int) [^2] | slivar_max_nhomalt | Maximum number of homozygous alternate alleles within population for small variants |  |
| String (Int) [^2] | slivar_max_ac | Maximum allele count within population for small variants |  |
| String (Int) [^2] | slivar_min_gq | Minimum genotype quality for small variants to be considered for compound heterozygous pairs |  |
| String | svpack_pop_vcfs | Comma-delimited array of structural variant population VCF paths |  |
| String | svpack_pop_vcf_indices | Comma-delimited array of structural variant population VCF index paths |  |

[^1]: Technically this value is interpreted as String by WDL, but slivar expects a Float, e.g, `0.03`.

[^2]: Technically these values are interpreted as String by WDL, but slivar expects an Int.
