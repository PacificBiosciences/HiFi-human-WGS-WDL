# bam_stats outputs

```wdl
{sample}.{movie}.read_length_and_quality.tsv.gz - per read length and quality metrics
```

## `{sample}.{movie}.read_length_and_quality.tsv.gz` - per read length and quality metrics

Base metrics are extracted for each read from the uBAM and stored in these 4 columns:

- movie
- read name
- read length: length of query sequence
- read quality: transformation of `rq` tag into Phred (log) space, e.g., `rq:f:0.99` (99% accuracy, 1 error in 100 bases) is Phred 20 ($-10 \times \log(1 - 0.99)$); this value is capped at Phred 60 for `rq:f:1.0`
