# family.wdl inputs and outputs

- [family.wdl inputs and outputs](#familywdl-inputs-and-outputs)
  - [DAG (simplified)](#dag-simplified)
  - [Inputs](#inputs)
    - [Family Struct](#family-struct)
    - [Sample Struct](#sample-struct)
  - [Outputs](#outputs)
    - [Alignments, Coverage, and QC](#alignments-coverage-and-qc)
    - [Small Variants (\<50 bp)](#small-variants-50-bp)
    - [Structural Variants (≥50 bp)](#structural-variants-50-bp)
    - [Tandem Repeat Genotyping](#tandem-repeat-genotyping)
    - [Variant Phasing](#variant-phasing)
    - [Variant Calling in Dark Regions](#variant-calling-in-dark-regions)
    - [5mCpG Methylation Calling](#5mcpg-methylation-calling)
    - [PGx Typing](#pgx-typing)
    - [Tertiary Analysis](#tertiary-analysis)

## DAG (simplified)

```mermaid
---
title: family.wdl
---
flowchart TD
  subgraph "create fail_reads bait FASTA"
    trgt_catalog["TRGT catalog BED"]
    bait_fasta["create bait FASTA"]
  end
  subgraph "`**Upstream of Phasing\n(per-sample)**`"
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
    samtools_merge["samtools merge"]
    mosdepth["mosdepth"]
    paraphase["Paraphase"]
    mitorsaw["MitorSaw"]
    deepvariant["DeepVariant"]
    sawfish_discover["Sawfish discover"]
  end
  subgraph "`**Joint Calling**`"
    glnexus["GLnexus (joint-call small variants)"]
    sawfish_call["Sawfish call"]
    split_glnexus["split small variant vcf by sample"]
    split_sawfish["split SV vcf by sample"]
  end
  subgraph "`**Phasing and Downstream**`"
    hiphase["HiPhase"]
    samtools_merge_fail_reads["samtools merge hifi_reads and fail_reads"]
    trgt["TRGT"]
    bam_stats["BAM stats"]
    bcftools_roh["bcftools roh"]
    bcftools_stats["bcftools stats\n(small variants)"]
    sv_stats["SV stats"]
    cpg_pileup["5mCpG pileup"]
    methbat["MethBat"]
    starphase["StarPhase"]
    pharmcat["PharmCat"]
  end
  subgraph " "
    merge_small_variants["bcftools merge small variants"]
    merge_svs["bcftools merge SV"]
    trgt_merge["trgt merge"]
  end
  subgraph "`**Tertiary Analysis**`"
    slivar_small_variants["slivar small variants"]
    svpack["svpack filter and annotate"]
    slivar_svpack["slivar svpack tsv"]
  end

  trgt_catalog --> bait_fasta --> bait_fail_reads
  fail_ubam --> bait_fail_reads --> pbmm2_align_fail_reads --> filter_fail_reads --> samtools_merge_fail_reads
  ubam --> pbmm2_align --> samtools_merge
  samtools_merge --> mosdepth
  samtools_merge --> paraphase
  samtools_merge --> mitorsaw
  samtools_merge_fail_reads --> trgt
  samtools_merge --> deepvariant
  samtools_merge --> sawfish_discover
  samtools_merge --> hiphase
  deepvariant --> sawfish_discover
  deepvariant --> glnexus
  sawfish_discover --> sawfish_call

  glnexus --> split_glnexus
  sawfish_call --> split_sawfish
  split_glnexus --> hiphase
  split_sawfish --> hiphase

  hiphase --> trgt
  hiphase --> bam_stats
  hiphase --> bcftools_roh
  hiphase --> bcftools_stats
  hiphase --> sv_stats
  hiphase --> cpg_pileup
  hiphase --> starphase
  hiphase --> pharmcat
  hiphase --> trgt_dropouts
  starphase --> pharmcat
  cpg_pileup --> methbat

  hiphase --> merge_small_variants
  hiphase --> merge_svs
  hiphase --> trgt_merge

  merge_small_variants --> slivar_small_variants
  merge_svs --> svpack
  svpack --> slivar_svpack
```

## Inputs

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| [Family](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/workflows/humanwgs_structs.wdl#L15) | family | Family struct describing samples, relationships, and unaligned BAM paths | [below](#family-struct) |
| File | [ref_map_file](./ref_map.md) | TSV containing reference genome file paths; must match backend | |
| String? | phenotypes | Comma-delimited list of HPO terms. | [Human Phenotype Ontology (HPO) phenotypes](https://hpo.jax.org/app/) associated with the cohort.<br/><br/>If omitted, tertiary analysis will be skipped. |
| File? | [tertiary_map_file](./tertiary_map.md) | TSV containing tertiary analysis file paths and thresholds; must match backend | `AF`/`AC`/`nhomalt` thresholds can be modified, but this will affect performance.<br/><br/>If omitted, tertiary analysis will be skipped. |
| Int | max_reads_per_alignment_chunk | Maximum reads per alignment chunk<br/><br/>Default: `500000` | |
| Int | pharmcat_min_coverage | Minimum coverage for PharmCAT<br/><br/>Default: `10` | |
| Int | glnexus_mem_gb | GLnexus memory<br/><br/>Default: `60` | |
| Boolean | use_gpu | Use GPU when possible<br/><br/>Default: `false` | [GPU support](./gpu.md#gpu-support) |
| Boolean | use_parabricks_deepvariant | Use Parabricks DeepVariant implementation<br/><br/>Default: `false` | If both `use_parabricks_deepvariant` and `use_gpu` are set to `true`, Parabricks DeepVariant will be used instead of standard DeepVariant.<br/><br/>[Parabricks DeepVariant](./parabricks.md#parabricks-deepvariant-subworkflow) |
| String | backend | Backend where the workflow will be executed<br/><br/>`["GCP", "Azure", "AWS-HealthOmics", "HPC"]` | |
| String? | zones | Zones where compute will take place; required if backend is set to 'GCP' | [Determining available zones in GCP](./backend-gcp#determining-available-zones) |
| String? | cpuPlatform | Minimum CPU platform to use for tasks on GCP | Optional, only necessary in certain zones lacking n1 nodes. |
| String? | gpuType | GPU type to use; required if use_gpu is set to `true` for cloud backends; must match backend | [Available GPU types](./gpu.md#gpu-types) |
| String? | container_registry | Container registry where workflow images are hosted.<br/><br/>Default: `"quay.io/pacbio"` | If omitted, [PacBio's public Quay.io registry](https://quay.io/organization/pacbio) will be used.<br/><br/>Custom container_registry must be set if backend is set to 'AWS-HealthOmics'. |
| Boolean | preemptible | Where possible, run tasks preemptibly<br/><br/>`[true, false]`<br/><br/>Default: `true` | If set to `true`, run tasks preemptibly where possible. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. |
| String? | debug_version | Debug version for testing purposes | |

### Family Struct

The `Family` struct contains the samples for the family. The struct has the following fields:

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | family_id | Unique identifier for the family | Alphanumeric characters, periods, dashes, and underscores are allowed. |
| Array\[[Sample](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/workflows/humanwgs_structs.wdl#L3)\] | samples | Sample struct with sample specific data and metadata. | [below](#sample-struct) |

### Sample Struct

The `Sample` struct contains sample specific data and metadata. The struct has the following fields:

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | sample_id | Unique identifier for the sample | Alphanumeric characters, periods, dashes, and underscores are allowed. |
| String? | sex | Sample sex<br/>`["MALE", "FEMALE", null]` | Used by HiFiCNV and TRGT for genotyping. Allosome karyotype will default to XX unless sex is specified as `"MALE"`.  Used for tertiary analysis X-linked inheritance filtering. |
| Boolean | affected | Affected status | If set to `true`, sample is described as being affected by all HPO terms in `phenotypes`.<br/>If set to `false`, sample is described as not being affected by all HPO terms in `phenotypes`. |
| Array\[File\] | hifi_reads | Array of paths to hifi_reads in unaligned BAM format. | |
| Array\[File\]? | fail_reads | Array of paths to fail_reads in unaligned BAM format (optional) | If provided, these reads will be aligned to the bait-captured regions. |
| String? | father_id | sample_id of father (optional) | |
| String? | mother_id | sample_id of mother (optional) | |

## Outputs

### Alignments, Coverage, and QC

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | workflow_name | Workflow name | |
| String | workflow_version | Workflow version | |
| Array\[String\] | msg | Messages from the workflow | |
| File | msg_file | File containing messages from the workflow | |
| Array\[String\] | sample_ids | Sample IDs | |
| File | stats_file | Table of summary statistics | |
| Array\[File\] | bam_statistics | BAM statistics | Per-read length and read-quality |
| Array\[File\] | read_length_plot | Distribution of read lengths | |
| Array\[File?\] | read_quality_plot | Distribution of read qualities | |
| Array\[File\] | merged_haplotagged_bam | Merged, haplotagged alignments | Includes unmapped reads |
| Array\[File\] | merged_haplotagged_bam_index | | |
| Array\[File\] | mosdepth_summary | Summary of aligned read depth | |
| Array\[File\] | mosdepth_region_bed | Median aligned read depth by 500bp windows | |
| Array\[File\] | mosdepth_region_bed_index | | |
| Array\[File\] | mosdepth_depth_distribution_plot | Distribution of aligned read depth | |
| Array\[File\] | mapq_distribution_plot | Distribution of mapping quality per alignment | |
| Array\[File\] | mg_distribution_plot | Distribution of gap-compressed identity per alignment | |
| Array\[String\] | stat_read_count | Number of reads | |
| Array\[String\] | stat_read_length_mean | Mean read length | |
| Array\[String\] | stat_read_length_median | Median read length | |
| Array\[String\] | stat_read_length_n50 | Read length N50 | |
| Array\[String\] | stat_read_quality_mean | Mean read quality | |
| Array\[String\] | stat_read_quality_median | Median read quality | |
| Array\[String\] | stat_mapped_read_count | Number of reads mapped to reference | |
| Array\[String\] | stat_mapped_read_percent | Percent of reads mapped to reference | |
| Array\[String\] | stat_gap_compressed_identity_mean | Mean gap-compressed identity | |
| Array\[String\] | stat_gap_compressed_identity_median | Median gap-compressed identity | |
| Array\[String\] | inferred_sex | Inferred sex | Sex is inferred based on relative depth of chrY alignments. |
| Array\[String\] | stat_depth_mean | Mean depth | |

### Small Variants (<50 bp)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File\] | phased_small_variant_vcf | Phased small variant VCF | |
| Array\[File\] | phased_small_variant_vcf_index | | |
| Array\[File?\] | small_variant_gvcf | Small variant GVCF | Can be used for joint-calling. |
| Array\[File?\] | small_variant_gvcf_index | | |
| Array\[File\] | small_variant_stats | Small variant statistics | Generated by `bcftools stats`. |
| Array\[String\] | stat_small_variant_SNV_count | Number of SNVs | (PASS variants) |
| Array\[String\] | stat_small_variant_INDEL_count | Number of INDELs | (PASS variants) |
| Array\[String\] | stat_small_variant_TSTV_ratio | Ts/Tv ratio | (PASS variants) |
| Array\[String\] | stat_small_variant_HETHOM_ratio | Het/Hom ratio for SNVs | (PASS variants) |
| Array\[File\] | snv_distribution_plot | Distribution of SNVs by REF, ALT | |
| Array\[File\] | indel_distribution_plot | Distribution of indels by size | |
| File? | joint_small_variants_vcf | Joint-called small variant VCF | |
| File? | joint_small_variants_vcf_index | | |

### Structural Variants (≥50 bp)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File\] | phased_sv_vcf | Phased structural variant VCF | |
| Array\[File\] | phased_sv_vcf_index | Index for phased structural variant VCF | |
| Array\[String\] | stat_sv_DUP_count | Number of DUP structural variants | (PASS variants) |
| Array\[String\] | stat_sv_DEL_count | Number of DEL structural variants | (PASS variants) |
| Array\[String\] | stat_sv_INS_count | Number of INS structural variants | (PASS variants) |
| Array\[String\] | stat_sv_INV_count | Number of INV structural variants | (PASS variants) |
| Array\[String\] | stat_sv_BND_count | Number of BND structural variants | (PASS variants) |
| Array\[String\] | stat_sv_SWAP_count | Number of structural variant sequence swap events | (PASS variants) |
| File | sv_supporting_reads | Supporting reads for structural variants | |
| Array\[File\] | sv_copynum_bedgraph | CNV copy number BEDGraph | |
| Array\[File\] | sv_depth_bw | CNV depth BigWig | |
| Array\[File\] | sv_gc_bias_corrected_depth_bw | CNV GC-bias corrected depth BigWig | |
| Array\[File\] | sv_maf_bw | CNV MAF BigWig | |
| Array\[File\] | sv_copynum_summary | CNV copy number summary JSON | |
| Array\[File\] | bcftools_roh_out | Regions of homozygosity | `bcftools roh` |
| Array\[File\] | bcftools_roh_bed | Regions of homozygosity BED | |
| File? | joint_sv_vcf | Joint-called structural variant VCF | |
| File? | joint_sv_vcf_index | | |

### Mitochondrial variants and haplotypes

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File\] | mitorsaw_vcf | Mitochondrial variant VCF | |
| Array\[File\] | mitorsaw_vcf_index | Index for mitochondrial variant VCF | |
| Array\[File\] | mitorsaw_hap_stats | Mitochondrial haplotype statistics | |

### Tandem Repeat Genotyping

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File\] | phased_trgt_vcf | Phased TRGT VCF | |
| Array\[File\] | phased_trgt_vcf_index | | |
| Array\[File\] | trgt_spanning_reads | Aligned TRGT spanning reads | |
| Array\[File\] | trgt_spanning_reads_index | | |
| Array\[File\] | trgt_coverage_dropouts | TRGT regions with coverage dropouts | |
| Array\[String\] | stat_trgt_genotyped_count | Number of sites genotyped by TRGT | |
| Array\[String\] | stat_trgt_uncalled_count | Number of sites ungenotyped by TRGT | |
| File? | joint_trgt_vcf | Joint-called TRGT VCF | |
| File? | joint_trgt_vcf_index | | |

### Variant Phasing

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File\] | phase_stats | Phasing statistics | |
| Array\[File\] | phase_blocks | Phase blocks | |
| Array\[File\] | phase_haplotags | Per-read phase assignment | |
| Array\[String\] | stat_phased_basepairs | Number of basepairs within phase blocks | |
| Array\[String\] | stat_phase_block_ng50 | Phase block NG50 | |

### Variant Calling in Dark Regions

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File?\] | paraphase_summary | Paraphase summary | |
| Array\[File?\] | paraphase_realigned_bam | BAM file of reads realigned by Paraphase | |
| Array\[File?\] | paraphase_realigned_bam_index | | |
| Array\[File?\] | paraphase_vcfs | Paraphase VCFs | Compressed as `.tar.gz` |

### 5mCpG Methylation Calling

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File?\] | cpg_hap1_bed | 5mCpG haplotype 1 BED | |
| Array\[File?\] | cpg_hap1_bed_index | | |
| Array\[File?\] | cpg_hap2_bed | 5mCpG haplotype 2 BED | |
| Array\[File?\] | cpg_hap2_bed_index | | |
| Array\[File?\] | cpg_combined_bed | 5mCpG combined BED | |
| Array\[File?\] | cpg_combined_bed_index | | |
| Array\[File?\] | cpg_hap1_bw | 5mCpG haplotype 1 BigWig | |
| Array\[File?\] | cpg_hap2_bw | 5mCpG haplotype 2 BigWig | |
| Array\[File?\] | cpg_combined_bw | 5mCpG combined BigWig | |
| Array\[String\] | stat_cpg_hap1_count | Number of scored reference 5mCpGs in haplotype 1 | |
| Array\[String\] | stat_cpg_hap2_count | Number of scored reference 5mCpGs in haplotype 2 | |
| Array\[String\] | stat_cpg_combined_count | Number of scored reference 5mCpGs combined | |
| Array\[File?\] | methbat_profile | MethBat 5mCpG profile | |
| Array\[String\] | stat_methbat_methylated_count | Number of profiled regions labeled as methylated | |
| Array\[String\] | stat_methbat_unmethylated_count | Number of profiled regions labeled as unmethylated | |
| Array\[String\] | stat_methbat_asm_count | Number of profiled regions labeled as having allele-specific methylation | |

### PGx Typing

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| Array\[File\] | pbstarphase_summary | StarPhase summary | Haplotype calls for PGx loci |
| Array\[File?\] | pharmcat_match_json | PharmCAT match JSON | |
| Array\[File?\] | pharmcat_phenotype_json | PharmCAT phenotype JSON | |
| Array\[File?\] | pharmcat_report_html | PharmCAT report HTML | |
| Array\[File?\] | pharmcat_report_json | PharmCAT report JSON | |

### Tertiary Analysis

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | tertiary_small_variant_filtered_vcf | Filtered, annotated small variant VCF | |
| File? | tertiary_small_variant_filtered_vcf_index | | |
| File? | tertiary_small_variant_filtered_tsv | Filtered, annotated small variant TSV | |
| File? | tertiary_small_variant_compound_het_vcf | Filtered, annotated compound heterozygous small variant VCF | |
| File? | tertiary_small_variant_compound_het_vcf_index | | |
| File? | tertiary_small_variant_compound_het_tsv | Filtered, annotated compound heterozygous small variant TSV | |
| File? | tertiary_sv_filtered_vcf | Filtered, annotated structural variant VCF | |
| File? | tertiary_sv_filtered_vcf_index | | |
| File? | tertiary_sv_filtered_tsv | Filtered, annotated structural variant TSV | |
