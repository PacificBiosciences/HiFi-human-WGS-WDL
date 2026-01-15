# TRGT repeat catalogs and inclusion of `fail_reads`

Our recommended TRGT Repeat Catalog is a merged set of ~70 pathogenic repeat sites ([STRchive](https://strchive.org)) and ~1M polymorphics repeats ([Adotto](https://zenodo.org/records/8329210)).  The file format is described in the [TRGT documentation](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md).

For some repeat loci, you may wish to include aligned `fail_reads` in addition to aligned `hifi_reads` for TRGT genotyping.  In this workflow, we mark these loci in the repeat catalog by adding the tag `INCLUDE_FAIL_READS` in column 4, as shown in the examples below.  If you provide `fail_reads`, we will align them to the reference genome and use them _**only**_ for TRGT genotyping, _**only**_ at these loci.  By default, we do not use `fail_reads` for TRGT genotyping, and we do not use `fail_reads` for any other purpose in this workflow.

Tagged loci in the default repeat catalog:

```text
chr4    39348424        39348483        ID=CANVAS_RFC1;MOTIFS=AAGGG,ACAGG,AGGGC,AAGGC,AGAGG,AAAAG,AAAGG,AAGAG,AAAGGG;STRUC=<TR>;INCLUDE_FAIL_READS
chr9    69037270        69037304        ID=FRDA_FXN;MOTIFS=A,GAA;STRUC=<TR>;INCLUDE_FAIL_READS
chr13   102161574       102161726       ID=SCA27B_FGF14;MOTIFS=GAA,GAAGGA,GAAGAAGAAGAAGCA,AAGGAG;STRUC=<TR>;INCLUDE_FAIL_READS
```

## Detecting coverage dropouts

Dropouts in coverage at TRGT catalog loci may indicate the presence of large expansions that are not fully spanned by HiFi reads.  To detect such dropouts, we run a script (`find_trgt_dropouts.py`) after TRGT genotyping to compares the observed coverage at each locus with a fixed threshold (2 reads per expected haplotype) and reports abnormal coverage in `*.trgt.dropouts.txt`.  The columns in this file are as follows:

- chrom
- start
- end
- trid
- expected_ploidy
- hap1_count
- hap2_count
- unphased_count
- fail_read_count
- dropout (FullDropout, HaplotypeDropout, PhasingDropout)

The `fail_read_count` column indicates the number of `fail_reads` that aligned to the locus.  This can help interpret dropouts at loci where `fail_reads` were included for genotyping.

The dropout column can be interpreted as follows:

- FullDropout: total HiFi read count (hap1 + hap2 + unphased) is less than `expected_ploidy * 2`
- HaplotypeDropout: one expected haplotype has fewer than 2 HiFi reads
- PhasingDropout: the total HiFi read depth is at least `expected_ploidy * 2`, but both haplotypes have fewer than 2 reads
