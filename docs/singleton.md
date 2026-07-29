# singleton.wdl inputs and outputs

- [singleton.wdl inputs and outputs](#singletonwdl-inputs-and-outputs)
  - [DAG (simplified)](#dag-simplified)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
    - [Alignments, Coverage, and QC](#alignments-coverage-and-qc)
    - [Small Variants (\<50 bp)](#small-variants-50-bp)
    - [Structural Variants (≥50 bp)](#structural-variants-50-bp)
    - [Tandem Repeat Genotyping](#tandem-repeat-genotyping)
    - [Variant Phasing](#variant-phasing)
    - [Variant Calling in Dark Regions](#variant-calling-in-dark-regions)
    - [Specialized Repeat Genotyping](#specialized-repeat-genotyping)
    - [5mCpG Methylation Calling](#5mcpg-methylation-calling)
    - [PGx Typing](#pgx-typing)

## DAG (simplified)

```mermaid
---
title: singleton.wdl
---
flowchart TD
  subgraph "create fail_reads bait FASTA"
    trgt_catalog["TRGT catalog BED"]
    bait_fasta["create bait FASTA"]
  end
  subgraph "`**Upstream of Phasing**`"
    subgraph "per hifi_reads uBAM"
      ubam[/"HiFi uBAM"/]
      pbmm2_align["pbmm2 align"]
    end
    subgraph "per fail_reads uBAM"
      fail_ubam[/"fail reads uBAM (if provided)"/]
      bait_fail_reads["baited fail reads (if fail_reads provided)"]
      pbmm2_align_fail_reads["pbmm2 align baited fail_reads (if fail_reads provided)"]
      filter_fail_reads["filter fail_reads alignments (if fail_reads provided)"]
    end
    pbsamoa_merge["pbsamoa merge"]
    mosdepth["mosdepth"]
    paraphase["Paraphase"]
    mitorsaw["MitorSaw"]
    kivvi_kiv2["kivvi (KIV2)"]
    kivvi_d4z4["kivvi (D4Z4)"]
    deepvariant["DeepVariant"]
    sawfish_discover["Sawfish discover"]
    sawfish_call["Sawfish call"]
  end
  subgraph "`**Phasing and Downstream**`"
    hiphase["HiPhase"]
    pbsamoa_merge_fail_reads["pbsamoa merge phased hifi_reads and aligned fail_reads"]
    trgt["TRGT"]
    pbjam["pbjam bam_stats"]
    bcftools_roh["bcftools roh"]
    bcftools_stats["bcftools stats\n(small variants)"]
    sv_stats["SV stats"]
    methbat_pileup["MethBat pileup"]
    methbat_profile["MethBat profile"]
    starphase["StarPhase"]
  end

  trgt_catalog --> bait_fasta --> bait_fail_reads
  fail_ubam --> bait_fail_reads --> pbmm2_align_fail_reads --> filter_fail_reads --> pbsamoa_merge_fail_reads
  ubam --> pbmm2_align --> pbsamoa_merge
  pbsamoa_merge --> mosdepth
  pbsamoa_merge --> paraphase
  pbsamoa_merge --> mitorsaw
  pbsamoa_merge --> kivvi_kiv2
  pbsamoa_merge --> kivvi_d4z4
  pbsamoa_merge_fail_reads --> trgt
  pbsamoa_merge --> deepvariant
  pbsamoa_merge --> sawfish_discover
  pbsamoa_merge --> hiphase
  deepvariant --> sawfish_discover
  deepvariant --> hiphase
  sawfish_discover --> sawfish_call --> hiphase

  hiphase --> trgt
  hiphase --> pbjam
  hiphase --> bcftools_roh
  hiphase --> bcftools_stats
  hiphase --> sv_stats
  hiphase --> methbat_pileup --> methbat_profile
  hiphase --> starphase
  hiphase --> trgt_dropouts
```

## Inputs

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | sample_id | Unique identifier for the sample | Alphanumeric characters, periods, dashes, and underscores are allowed. |
| Array\[File\] | hifi_reads | Array of paths to hifi_reads in unaligned BAM format. | |
| Array\[File\]? | fail_reads | Array of paths to fail_reads in unaligned BAM format (optional) | If provided, these reads will be aligned to the bait-captured regions. |
| String | ref_name | Reference genome to use for this workflow run<br/><br/>`["GRCh38", "GRCh38_GIABv3"]`<br/><br/>Default: `"GRCh38"` | |
| File? | trgt_tandem_repeat_bed_override | Optional BED file to override the default TRGT tandem repeat catalog | |
| File? | methbat_region_tsv_override | Optional TSV file to override the default MethBat methylation profiling regions | |
| Boolean | use_alignment_chunking | Whether to chunk BAM files for alignment. If false, all reads will be aligned in a single chunk.<br/><br/>Default: `true` | |
| Boolean | use_gpu | Use GPU when possible<br/><br/>Default: `false` | [GPU support](./gpu.md#gpu-support) |
| Boolean | use_parabricks_deepvariant | Use Parabricks DeepVariant implementation<br/><br/>Default: `false` | If both `use_parabricks_deepvariant` and `use_gpu` are set to `true`, Parabricks DeepVariant will be used instead of standard DeepVariant.<br/><br/>[Parabricks DeepVariant](./parabricks.md#parabricks-deepvariant-subworkflow) |
| String | backend | Backend where the workflow will be executed<br/><br/>`["GCP", "Azure", "AWS-HealthOmics", "HPC"]` | |
| String? | zones | Zones where compute will take place; required if backend is set to 'GCP' | [Determining available zones in GCP](./backend-gcp.md#determining-available-zones) |
| String? | cpuPlatform | Minimum CPU platform to use for tasks on GCP | Optional, only necessary in certain zones lacking n1 nodes. |
| String? | gpuType | GPU type to use; required if use_gpu is set to `true` for cloud backends; must match backend | [Available GPU types](./gpu.md#gpu-types) |
| String? | container_registry | Container registry where workflow images are hosted.<br/><br/>Default: `"quay.io/pacbio"` | If omitted, [PacBio's public Quay.io registry](https://quay.io/organization/pacbio) will be used.<br/><br/>Custom container_registry must be set if backend is set to 'AWS-HealthOmics'. |
| Boolean | preemptible | Where possible, run tasks preemptibly<br/><br/>`[true, false]`<br/><br/>Default: `true` | If set to `true`, run tasks preemptibly where possible. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. |
| String? | debug_version | Debug version for testing purposes | |

## Outputs

### Alignments, Coverage, and QC

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | workflow_name | Workflow name | |
| String | workflow_version | Workflow version | |
| Array\[String\] | msg | Messages from the workflow | |
| File | msg_file | File containing messages from the workflow | |
| File | stats_file | Table of summary statistics | |
| File | read_length_plot | Distribution of read lengths | |
| File? | read_quality_plot | Distribution of read qualities | |
| File | merged_haplotagged_bam | Merged, haplotagged alignments | Includes unmapped reads |
| File | merged_haplotagged_bam_index | | |
| File | mosdepth_summary | Summary of aligned read depth | |
| File | mosdepth_region_bed | Median aligned read depth by 500bp windows | |
| File | mosdepth_region_bed_index | | |
| File | mosdepth_depth_distribution_plot | Distribution of aligned read depth | |
| File | mapq_distribution_plot | Distribution of mapping quality per alignment | |
| File | mg_distribution_plot | Distribution of gap-compressed identity per alignment | |
| String | stat_read_count | Number of reads | |
| String | stat_read_length_mean | Mean read length | |
| String | stat_read_length_median | Median read length | |
| String | stat_read_length_n50 | Read length N50 | |
| String | stat_read_quality_mean | Mean read quality | |
| String | stat_read_quality_median | Median read quality | |
| String | stat_mapped_read_count | Number of reads mapped to reference | |
| String | stat_mapped_read_percent | Percent of reads mapped to reference | |
| String | stat_gap_compressed_identity_mean | Mean gap-compressed identity | |
| String | stat_gap_compressed_identity_median | Median gap-compressed identity | |
| String | inferred_sex | Inferred sex | Sex is inferred based on relative depth of chrY alignments. |
| String | stat_depth_mean | Mean depth | |

### Small Variants (<50 bp)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phased_small_variant_vcf | Phased small variant VCF | |
| File | phased_small_variant_vcf_index | | |
| File? | small_variant_gvcf | Small variant GVCF | Can be used for joint-calling. |
| File? | small_variant_gvcf_index | | |
| File | small_variant_stats | Small variant statistics | Generated by `bcftools stats`. |
| String | stat_small_variant_SNV_count | Number of SNVs | (PASS variants) |
| String | stat_small_variant_INDEL_count | Number of INDELs | (PASS variants) |
| String | stat_small_variant_TSTV_ratio | Ts/Tv ratio | (PASS variants) |
| String | stat_small_variant_HETHOM_ratio | Het/Hom ratio for SNVs | (PASS variants) |
| File | snv_distribution_plot | Distribution of SNVs by REF, ALT | |
| File | indel_distribution_plot | Distribution of indels by size | |

### Structural Variants (≥50 bp)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phased_sv_vcf | Phased structural variant VCF | |
| File | phased_sv_vcf_index | Index for phased structural variant VCF | |
| String | stat_sv_DUP_count | Number of DUP structural variants | (PASS variants) |
| String | stat_sv_DEL_count | Number of DEL structural variants | (PASS variants) |
| String | stat_sv_INS_count | Number of INS structural variants | (PASS variants) |
| String | stat_sv_INV_count | Number of INV structural variants | (PASS variants) |
| String | stat_sv_BND_count | Number of BND structural variants | (PASS variants) |
| String | stat_sv_SWAP_count | Number of structural variant sequence swap events | (PASS variants) |
| File | sv_stats_plot | Distribution of DEL/INS/DUP/INV by size | |
| File | sv_supporting_reads | Supporting reads for structural variants | |
| File | sv_copynum_bedgraph | CNV copy number BEDGraph | |
| File | sv_depth_bw | CNV depth BigWig | |
| File | sv_gc_bias_corrected_depth_bw | CNV GC-bias corrected depth BigWig | |
| File | sv_maf_bw | CNV MAF BigWig | |
| File | sv_copynum_summary | CNV copy number summary JSON | |
| File | bcftools_roh_out | Regions of homozygosity | `bcftools roh` |
| File | bcftools_roh_bed | Regions of homozygosity BED | |

### Mitochondrial variants and haplotypes

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | mitorsaw_vcf | Mitochondrial variant VCF | |
| File | mitorsaw_vcf_index | Index for mitochondrial variant VCF | |
| File | mitorsaw_hap_stats | Mitochondrial haplotype statistics | |

### Tandem Repeat Genotyping

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phased_trgt_vcf | Phased TRGT VCF | |
| File | phased_trgt_vcf_index | | |
| File | trgt_spanning_reads | Aligned TRGT spanning reads | |
| File | trgt_spanning_reads_index | | |
| File | trgt_coverage_dropouts | TRGT regions with coverage dropouts | |
| String | stat_trgt_genotyped_count | Number of sites genotyped by TRGT | |
| String | stat_trgt_uncalled_count | Number of sites ungenotyped by TRGT | |

### Variant Phasing

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phase_stats | Phasing statistics | |
| File | phase_blocks | Phase blocks | |
| File | phase_haplotags | Per-read phase assignment | |
| String | stat_phased_basepairs | Number of basepairs within phase blocks | |
| String | stat_phase_block_ng50 | Phase block NG50 | |

### Variant Calling in Dark Regions

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | paraphase_summary | Paraphase summary | |
| File? | paraphase_realigned_bam | BAM file of reads realigned by Paraphase | |
| File? | paraphase_realigned_bam_index | | |
| File? | paraphase_vcfs | Paraphase VCFs | Compressed as `.tar.gz` |

### Specialized Repeat Genotyping

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | kivvi_kiv2_vcf | KIV2 repeat variant VCF | |
| File? | kivvi_kiv2_vcf_index | | |
| File? | kivvi_kiv2_json | KIV2 repeat genotype JSON | |
| File? | kivvi_kiv2_realigned_bam | KIV2-realigned BAM | |
| File? | kivvi_kiv2_realigned_bam_index | | |
| File? | kivvi_kiv2_allele_plot | KIV2 assembled allele plot | |
| File? | kivvi_d4z4_vcf | D4Z4 repeat variant VCF | |
| File? | kivvi_d4z4_vcf_index | | |
| File? | kivvi_d4z4_json | D4Z4 repeat genotype JSON | |
| File? | kivvi_d4z4_realigned_bam | D4Z4-realigned BAM | |
| File? | kivvi_d4z4_realigned_bam_index | | |
| File? | kivvi_d4z4_allele_plot | D4Z4 assembled allele plot | |

### 5mCpG Methylation Calling

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | cpg_pileup_bed | 5mCpG pileup BED | |
| File? | cpg_pileup_bed_index | | |
| File? | hmcpg_pileup_bed | 5hmCpG pileup BED | |
| File? | hmcpg_pileup_bed_index | | |
| String | stat_cpg_hap1_count | Number of scored reference 5mCpGs in haplotype 1 | |
| String | stat_cpg_hap2_count | Number of scored reference 5mCpGs in haplotype 2 | |
| String | stat_cpg_combined_count | Number of scored reference 5mCpGs combined | |
| File? | methbat_profile | MethBat 5mCpG profile | |
| String | stat_methbat_methylated_count | Number of profiled regions labeled as methylated | |
| String | stat_methbat_unmethylated_count | Number of profiled regions labeled as unmethylated | |
| String | stat_methbat_asm_count | Number of profiled regions labeled as having allele-specific methylation | |

### PGx Typing

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | pbstarphase_summary | StarPhase summary | Haplotype calls for PGx loci |
| File? | pbstarphase_tsv | StarPhase summary in TSV format for PharmCAT | |
