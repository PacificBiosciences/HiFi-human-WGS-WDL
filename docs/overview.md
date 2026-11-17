# Pipeline overview

This describes what the workflow does end to end, stage by stage, and why. For the full input/output list and a generated call-graph DAG, see [docs/singleton.md](singleton.md) and [docs/family.md](family.md). For the full tool list, citations, and pinned versions, see [Key Tools](tools.md) and [Tool versions and Containers](tools_containers.md).

There are two entrypoints:

- **`singleton.wdl`** analyzes one sample.
- **`family.wdl`** analyzes a family of related samples (a `Family` struct: a `family_id` plus an array of `Sample`s, each a `sample_id` with its `hifi_reads`/`fail_reads`). It runs every sample through the same per-sample stages as `singleton.wdl`, adds a joint-calling stage across the family, then merges family-wide VCFs at the end.

`family.wdl` is a superset of `singleton.wdl`: everything below through [phasing and per-sample outputs](#5-phasing-hiphase) applies to both. The [joint calling](#joint-calling-familywdl-only) and [family-wide merge](#family-wide-merge-familywdl-only) sections apply only to `family.wdl`, and only when the family has more than one sample (`single_sample` is derived from `length(family.samples) == 1` and short-circuits both).

## 1. Reference data and TRGT catalog setup

Before any sample is touched, the workflow unpacks a pinned [reference data container](reference_data_containers.md) ([doi:10.5281/zenodo.21517827](https://doi.org/10.5281/zenodo.21517827)) (keyed by the `ref_name` input, e.g. `GRCh38` or `GRCh38_GIABv3`) into a reference FASTA plus the region files every downstream tool needs: the default TRGT tandem-repeat catalog, Sawfish's exclude/expected-copy-number BEDs, and the MethBat profiling regions. `trgt_tandem_repeat_bed_override`/`methbat_region_tsv_override` inputs replace the packaged defaults if supplied.

The TRGT catalog is then scanned once for any regions flagged `INCLUDE_FAIL_READS`, and a bait FASTA is built from them -- used in step 2 to selectively rescue `fail_reads` that land in those regions, rather than aligning every failed read genome-wide.

## 2. Alignment

Each element of the `hifi_reads` array (unaligned HiFi BAM) is aligned independently with **[pbmm2](pbmm2.md)** ([PacificBiosciences/pbmm2](https://github.com/PacificBiosciences/pbmm2), PacBio's `minimap2`-based aligner tuned for HiFi reads; `use_alignment_chunking` controls whether an individual BAM is split into chunks for parallel alignment). If more than one `hifi_reads` element was supplied, the resulting alignments are merged into one BAM with **[pbsamoa](https://github.com/PacificBiosciences/pbsamoa)** (a C++ library for SAM/BAM/CRAM I/O).

If `fail_reads` (reads that failed HiFi quality filtering) are also supplied, they're separately aligned against just the bait-captured regions from step 1, filtered to keep only reads that actually land in an `INCLUDE_FAIL_READS` region, and merged if there's more than one element. This is a targeted rescue, not a full realignment of `fail_reads` -- it only recovers reads useful for tandem-repeat genotyping (step 6) in specific dropout-prone regions, and is skipped entirely if `fail_reads` isn't provided.

### Optional: nucleosome and MSP tagging (`run_fibertools`)

If `run_fibertools` is `true` (default `false`), **[fibertools](https://fiberseq.github.io/)** ([Jha et al., 2024](https://doi.org/10.1101/gr.279095.124)) `add-nucleosomes` runs on the merged `hifi_reads` alignment and adds nucleosome and methyltransferase-sensitive patch (MSP) annotations to each read, using the m6A base-modification predictions already in the BAM. It is intended for Fiber-seq data; the input `hifi_reads` must carry m6A modification tags.

The nucleosome-tagged BAM is a separate file; the plain aligned BAM is kept. The tagged BAM replaces the plain one for [phasing](#5-phasing-hiphase) and everything after it (TRGT, MethBat, StarPhase, QC stats). Steps 3 and 4 (mosdepth, Paraphase, MitorSaw, Kivvi, DeepVariant/Parabricks, Sawfish `discover` and joint `call`) always read the untagged BAM.

## 3. Coverage, sex inference, and per-sample dark-region/repeat calling

Once a sample's HiFi reads are aligned, several independent tasks run off that one alignment:

- **[mosdepth](https://github.com/brentp/mosdepth)** ([Pedersen & Quinlan, 2018](https://doi.org/10.1093/bioinformatics/btx699)) computes depth summaries and infers sex from relative chrX/chrY depth (compared against `max_norm_female_chrY_depth` from the reference bundle); this inferred value is what every downstream step (TRGT, Sawfish) uses.
- **[Paraphase](https://github.com/PacificBiosciences/paraphase)** ([Chen et al., 2025](https://doi.org/10.1038/s41467-025-57505-2)) calls variants in paralogous gene families with high sequence similarity (e.g. _SMN1_/_SMN2_, _PMS2_/_PMS2CL_) that standard alignment-based calling handles poorly.
- **[MitorSaw](https://github.com/PacificBiosciences/mitorsaw)** calls mitochondrial variants and haplotype statistics from the (nuclear) alignment.
- **[Kivvi](https://github.com/PacificBiosciences/kivvi)** ([Chen et al., 2026, preprint](https://doi.org/10.64898/2026.04.10.717730)) genotypes two specific hard-to-call repeat loci directly: KIV2 (LPA, cardiovascular risk) and D4Z4 (FSHD).

None of these four block or feed each other -- they all read the same aligned BAM independently.

## 4. Variant calling

- **Small variants (<50 bp):** **[DeepVariant](deepvariant.md)** ([documentation](https://github.com/google/deepvariant/tree/master/docs); [Poplin et al., 2018](https://doi.org/10.1038/nbt.4235)) calls SNVs/indels from the aligned BAM, running its three internal steps (`make_examples` / `call_variants` / `postprocess_variants`) separately so each can be sized for CPU or GPU independently. If `use_gpu` and `use_parabricks_deepvariant` are both true, this is replaced wholesale by **[Parabricks DeepVariant](parabricks.md)** ([NVIDIA docs](https://docs.nvidia.com/clara/parabricks/latest/)) instead. See [GPU support](gpu.md). DeepVariant calls can differ slightly between hardware types; see [How deterministic are variant calling results?](determinism.md).
- **Structural variants (≥50 bp):** **[Sawfish](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md)** ([Saunders et al., 2025](https://doi.org/10.1093/bioinformatics/btaf136)) `discover` builds a per-sample structural-variant evidence tarball from the aligned BAM. What happens next depends on the entrypoint:
  - `singleton.wdl` (and any single-sample `family.wdl` run): Sawfish `call` genotypes directly from that one sample's `discover` tarball, right here in the per-sample stage.
  - `family.wdl` with more than one sample: genotyping is deferred to the family-wide [joint calling](#joint-calling-familywdl-only) stage below, so SVs can be called with evidence from every family member at once.

## Joint calling (`family.wdl` only)

Before phasing, when a family has more than one sample, small variants and SVs are called jointly across the whole family instead of per sample, then split back into per-sample VCFs for phasing:

- **[GLnexus](https://github.com/dnanexus-rnd/GLnexus/wiki)** ([Yun et al., 2020](https://doi.org/10.1093/bioinformatics/btaa1081)) joint-genotypes small variants from every sample's DeepVariant gVCF into one family-level callset.
- Sawfish `call` genotypes SVs from every sample's `discover` tarball and aligned BAM together, producing family-wide copy-number tracks (BEDGraph/BigWig) alongside the VCF.

Both joint VCFs are then split by sample ([`bcftools`](http://www.htslib.org/doc/bcftools.html), [Danecek et al., 2021](https://doi.org/10.1093/gigascience/giab008)) before each sample proceeds to phasing (step 5) independently. Joint calling is why `family.wdl` finds SVs and small variants that per-sample calling alone would miss or genotype less confidently -- evidence from relatives informs each sample's genotype -- at the cost of needing every family member's data before any sample can be phased.

## 5. Phasing (HiPhase)

**[HiPhase](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md)** ([Holt et al., 2024](https://doi.org/10.1093/bioinformatics/btae042)) phases the small-variant and structural-variant VCFs for a sample together against its own aligned HiFi reads, in one joint step (small variants and SVs share phase information, which produces more consistent phase blocks than phasing them separately). It also haplotags the alignment itself, producing the `merged_haplotagged_bam` that most downstream steps and the final output use in place of the plain aligned BAM. If `run_fibertools` is enabled, the BAM HiPhase haplotags is the nucleosome-tagged one, so the haplotagged BAM carries the nucleosome/MSP annotations as well.

For `family.wdl`, this step runs per sample, using that sample's split-out joint-called VCFs (see above) rather than its own unmerged per-sample calls.

## 6. Tandem repeat genotyping (TRGT)

**[TRGT](https://github.com/PacificBiosciences/trgt)** ([Dolzhenko et al., 2024](https://doi.org/10.1038/s41587-023-02057-3)) genotypes the tandem-repeat catalog (step 1) against the haplotagged, phased alignment, plus any rescued `fail_reads` from step 2 as extra spanning-read support in dropout-prone regions. Output is phased where phase information from step 5 is available for a given locus.

## 7. Methylation profiling (MethBat)

If the input `hifi_reads` carry base-modification (`MM`/`ML`) tags, **[MethBat](https://github.com/PacificBiosciences/MethBat)** produces a 5mCpG (and 5hmCpG, if present) pileup from the haplotagged alignment, then profiles it against the MethBat region set from step 1 to call methylation status (methylated / unmethylated / allele-specific) per region. If no modification tags are present, these outputs are simply omitted -- not an error.

## 8. HLA and PGx typing (StarPhase)

**[StarPhase](https://github.com/PacificBiosciences/pb-StarPhase/blob/main/docs/user_guide.md)** ([Holt et al., 2024, preprint](https://doi.org/10.1101/2024.12.10.627527)) calls HLA and pharmacogenomic (PGx) diplotypes from the phased small-variant VCF and alignment. Whether it runs is controlled per-reference by `run_starphase`, baked into the reference data container from step 1 -- in practice true for both currently supported GRCh38-based `ref_name` options (`GRCh38`, `GRCh38_GIABv3`), so this step runs by default.

## 9. QC stats and reporting

`bcftools stats`/`roh` ([Narasimhan et al., 2016](https://doi.org/10.1093/bioinformatics/btw044), for the `roh` algorithm specifically), [`pbjam`](https://github.com/PacificBiosciences/pbjam) `bam_stats`, and per-category SV-size stats all summarize the outputs above (read/alignment QC, variant counts, regions of homozygosity, SV size distribution) into the plots and summary tables listed in [docs/singleton.md](singleton.md#outputs)/[docs/family.md](family.md#outputs). Warnings raised by individual tasks along the way (e.g. alignment, TRGT catalog processing, Kivvi, Sawfish discover, TRGT genotyping) are collected into the `msg`/`msg_file` outputs.

## Family-wide merge (`family.wdl` only)

When a family has more than one sample, after every sample has been phased independently (step 5), their phased small-variant VCFs, SV VCFs, and TRGT VCFs are each merged (`bcftools merge` for variants, `trgt merge` for repeats) into one family-wide VCF per variant type, in addition to the per-sample outputs each sample already produced on its own.

## Singleton vs. family, summarized

| | `singleton.wdl` | `family.wdl` |
| --- | --- | --- |
| Samples per run | 1 | 1 or more, related |
| Small/structural variant calling | Per sample only | Per sample if 1 sample; joint (GLnexus/Sawfish) across the family otherwise |
| Family-wide merged VCFs | N/A | Produced when there's more than 1 sample |
| Per-sample outputs (alignment, TRGT, methylation, PGx, mito, dark-region/repeat calling) | Yes | Yes, for every sample |
