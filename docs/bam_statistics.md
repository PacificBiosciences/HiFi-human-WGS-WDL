# bam_stats outputs

## `bam_statistics`

A compressed TSV file with a row for each record in the haplotagged BAM and the following columns:

- movie name
- read name
- read length
- Phred scaled read quality
- alignment type (unmapped, primary, supplementary; because supplementary alignments are included, reads may appear on multiple rows)
- mapping quality (`MAPQ`), if mapped
- gap-compressed identity (`mg`), if mapped

## `read_length_plot`

A histogram of read lengths, using only records marked `prim` or `unmapped`.

## `read_quality_plot`

A histogram of read qualities, using only records marked `prim` or `unmapped`. This output is only generated if the input BAMs contain the `rq` tag.

## `mapq_distribution_plot`, `mg_distribution_plot`

A histogram of mapping qualities and gap-compressed identities, respectively.

## `stat_num_reads`, `stat_read_length_mean`, `stat_read_length_median`, `stat_read_length_n50`, `stat_read_quality_mean`, `stat_read_quality_median`

Statistics computed using only records marked `prim` or `unmapped`.

## `stat_mapped_read_count`, `stat_mapped_percent`

Count of primary alignments, and primary alignments as a percentage of total reads.

## `stat_mean_gap_compressed_identity`

Mean gap-compressed identity of primary and supplementary alignments.
