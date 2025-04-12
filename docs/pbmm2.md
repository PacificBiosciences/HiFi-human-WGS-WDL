# pbmm2 alignment Subworkflow

```mermaid
flowchart TD
  hifi_reads[/"HiFi reads BAm"/] --> is_aligned{"is aligned?"}
  is_aligned -- yes --> samtools_reset["samtools reset"]
  is_aligned -- no --> has_kinetics{"kinetics?"}
  has_kinetics -- yes --> samtools_reset
  has_kinetics -- no --> count_records["count records"]
  samtools_reset --> count_records
  count_records --> compare_counts{"compare counts?"}
  compare_counts -- yes --> chunk_bam["chunk BAM"]
  compare_counts -- no --> pbmm2_align["pbmm2 align"]
  chunk_bam --> pbmm2_align
```

This subworkflow checks an input BAM for evidence of alignment or kinetics. If it finds either of these, it strips alignment and kinetics information.  Next, it counts the number of records in the BAM, and if  chunking is enabled and the number of records is greater than `max_reads_per_chunk`, the BAM is split into chunks of no larger than `max_reads_per_chunk`. Finally, chunks are aligned to the reference with pbmm2.
