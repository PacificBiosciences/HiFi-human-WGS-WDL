# bam_stats outputs

```wdl
{sample}.{movie}.read_length_and_quality.tsv - per read length and quality metrics
{sample}.{movie}.read_length_summary.tsv - histogram of read lengths
{sample}.{movie}.read_quality_summary.tsv - histogram of read qualities
```

## `{sample}.{movie}.read_length_and_quality.tsv.gz` - per read length and quality metrics

Base metrics are extracted for each read from the uBAM and stored in these 4 columns:

- movie
- read name
- read length: length of query sequence
- read quality: transformation of `rq` tag into Phred (log) space, e.g., `rq:f:0.99` (99% accuracy, 1 error in 100 bases) is Phred 20 ($-10 \times \log(1 - 0.99)$); this value is capped at Phred 60 for `rq:f:1.0`

## `{sample}.{movie}.*_summary.tsv` - read length and quality histograms

Values are binned for the histogram files.

For a given row, the values are:

- column 1: bin name; lower limit of bin; reads counted in this bin fall into the interval `[row N column 1, row N+1 column 1)`
- column 2: count of reads in this bin
- column 3: count of base pairs in this bin

As an example these are the first few rows of a read_length_summary.tsv:

```tsv
0       38      28079
1000    284     466848
2000    925     2383460
```

First row are reads with 0-999 bp/read. There are 38 reads in this bin. They sum to 28 kbp.
Second row are reads with 1000-1999 bp/read. There are 284 reads in this bin. They sum to 466 kbp.
Third row are reads with 2000-2999 bp/read. There are 925 reads in this bin. They sum to 2.38 Mbp.
And here are some rows from the middle of a read_quality_summary.tsv:

```tsv
28      370375  6973927935
29      387446  7256330089
30      372960  6888342745
```

First row are reads with Phred scaled read quality in `[28, 29)`. These have predicted error rates between ~1/631 (Phred 28, inclusive) and ~1/794 (Phred 29). There are 370k reads in this bin, and they sum to 6.97 Gbp.
