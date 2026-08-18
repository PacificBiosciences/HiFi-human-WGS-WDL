# Key Tools

This page lists the primary bioinformatics tools this pipeline relies on for
its scientific analysis, with a link to each tool's own documentation. For
every pinned tool/library version (including supporting libraries like
samtools and Python packages) and container digests, see [Tool versions and
Containers](tools_containers.md).

## Alignment

- **pbmm2** -- PacBio's `minimap2`-based aligner tuned for HiFi reads.
  [Documentation](https://github.com/PacificBiosciences/pbmm2)
- **Minimap2** -- general-purpose long-read aligner; used internally by
  [Paraphase](#paralog-genes) for realignment.
  [Documentation](https://lh3.github.io/minimap2/)

## Coverage

- **mosdepth** -- fast BAM/CRAM depth calculation; used for depth summaries
  and sex inference.
  [Documentation](https://github.com/brentp/mosdepth)

## Small variants

- **DeepVariant** -- deep-learning small variant caller (SNVs/indels).
  [Documentation](https://github.com/google/deepvariant/tree/master/docs)
- **NVIDIA Parabricks** -- GPU-accelerated DeepVariant pipeline, used when
  `use_parabricks_deepvariant=true`.
  [Documentation](https://docs.nvidia.com/clara/parabricks/latest/)
- **GLnexus** -- joint genotyper that merges per-sample GVCFs into a
  family-level callset (`family.wdl` only).
  [Documentation](https://github.com/dnanexus-rnd/GLnexus/wiki)

## Structural variants

- **Sawfish** -- structural variant and copy-number caller.
  [Documentation](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md)

## Tandem repeats

- **TRGT** -- targeted tandem repeat genotyper.
  [Documentation](https://github.com/PacificBiosciences/trgt)

## Phasing

- **HiPhase** -- haplotype phasing across small variants, structural
  variants, and tandem repeats.
  [Documentation](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md)

## Methylation

- **MethBat** -- 5mC/5hmC methylation pileup and region-level profiling.
  [Documentation](https://github.com/PacificBiosciences/MethBat) (see also its [profile guide](https://github.com/PacificBiosciences/MethBat/blob/main/docs/profile_guide.md))

## Mitochondrial variants

- **MitorSaw** -- mitochondrial variant caller and haplotype statistics.
  [Documentation](https://github.com/PacificBiosciences/mitorsaw)

## Paralog genes

- **Paraphase** -- paralog-aware variant calling for genes with high
  sequence similarity (e.g., _SMN1_/_SMN2_, _PMS2_/_PMS2CL_).
  [Documentation](https://github.com/PacificBiosciences/paraphase)
- **Kivvi** -- targeted genotyping for the KIV2 (LPA) and D4Z4 (FSHD) repeat
  regions.
  [Documentation](https://github.com/PacificBiosciences/kivvi)

## Pharmacogenomics

- **StarPhase** -- HLA typing and pharmacogenomic (PGx) diplotype calling.
  [Documentation](https://github.com/PacificBiosciences/pb-StarPhase/blob/main/docs/user_guide.md)
