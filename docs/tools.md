# Key Tools

This page lists the primary bioinformatics tools this pipeline relies on for its scientific analysis, with a link to each tool's own documentation and a citation -- the tool's own publication where one exists, otherwise its git repository. For every pinned tool/library version and container digests, see [Tool versions and Containers](tools_containers.md).

## Alignment

- **pbmm2** -- PacBio's `minimap2`-based aligner tuned for HiFi reads.<br>
  [Documentation](https://github.com/PacificBiosciences/pbmm2)<br>
  Citation: [PacificBiosciences. pbmm2. GitHub repository.](https://github.com/PacificBiosciences/pbmm2)
- **Minimap2** -- general-purpose long-read aligner; used internally by [Paraphase](#paralog-genes) for realignment.<br>
  [Documentation](https://lh3.github.io/minimap2/)<br>
  Citation: [Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100.](https://doi.org/10.1093/bioinformatics/bty191)
- **pbsamoa** -- C++ library for SAM/BAM/CRAM I/O; used here to merge aligned BAMs (multi-chunk `hifi_reads` alignments, rescued `fail_reads` alignments).<br>
  [Documentation](https://github.com/PacificBiosciences/pbsamoa)<br>
  Citation: [PacificBiosciences. pbsamoa. GitHub repository.](https://github.com/PacificBiosciences/pbsamoa)

## Coverage

- **mosdepth** -- fast BAM/CRAM depth calculation; used for depth summaries and sex inference.<br>
  [Documentation](https://github.com/brentp/mosdepth)<br>
  Citation: [Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867-868.](https://doi.org/10.1093/bioinformatics/btx699)

## Small variants

- **DeepVariant** -- deep-learning small variant caller (SNVs/indels).<br>
  [Documentation](https://github.com/google/deepvariant/tree/master/docs)<br>
  Citation: [Poplin, R., et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36, 983-987.](https://doi.org/10.1038/nbt.4235)
- **NVIDIA Parabricks** -- GPU-accelerated DeepVariant pipeline, used when `use_parabricks_deepvariant=true`.<br>
  [Documentation](https://docs.nvidia.com/clara/parabricks/latest/)<br>
  Citation: [NVIDIA. Clara Parabricks documentation.](https://docs.nvidia.com/clara/parabricks/latest/) (no public git repository; also see the DeepVariant citation above)
- **GLnexus** -- joint genotyper that merges per-sample GVCFs into a family-level callset (`family.wdl` only).<br>
  [Documentation](https://github.com/dnanexus-rnd/GLnexus/wiki)<br>
  Citation: [Yun, T., et al. (2020). Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Bioinformatics, 36(24), 5582-5589.](https://doi.org/10.1093/bioinformatics/btaa1081)

## Structural variants

- **Sawfish** -- structural variant and copy-number caller.<br>
  [Documentation](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md)<br>
  Citation: [Saunders, C. T., et al. (2025). Sawfish: improving long-read structural variant discovery and genotyping with local haplotype modeling. Bioinformatics.](https://doi.org/10.1093/bioinformatics/btaf136)

## Tandem repeats

- **TRGT** -- targeted tandem repeat genotyper.<br>
  [Documentation](https://github.com/PacificBiosciences/trgt)<br>
  Citation: [Dolzhenko, E., et al. (2024). Characterization and visualization of tandem repeats at genome scale. Nature Biotechnology, 42, 1606-1614.](https://doi.org/10.1038/s41587-023-02057-3)

## Phasing

- **HiPhase** -- haplotype phasing across small variants, structural variants, and tandem repeats.<br>
  [Documentation](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md)<br>
  Citation: [Holt, J. M., et al. (2024). HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing. Bioinformatics, 40(2), btae042.](https://doi.org/10.1093/bioinformatics/btae042)

## Chromatin (Fiber-seq)

- **fibertools** -- adds nucleosome and methyltransferase-sensitive patch (MSP) annotations to reads carrying m6A predictions (`ft add-nucleosomes`); runs only when `run_fibertools=true`.<br>
  [Documentation](https://fiberseq.github.io/)<br>
  Citation: [Jha, A., et al. (2024). DNA-m6A calling and integrated long-read epigenetic and genetic analysis with fibertools. Genome Research.](https://doi.org/10.1101/gr.279095.124)

## Methylation

- **MethBat** -- 5mC/5hmC methylation pileup and region-level profiling.<br>
  [Documentation](https://github.com/PacificBiosciences/MethBat) (see also its [profile guide](https://github.com/PacificBiosciences/MethBat/blob/main/docs/profile_guide.md))<br>
  Citation: [PacificBiosciences. MethBat. GitHub repository.](https://github.com/PacificBiosciences/MethBat)

## Mitochondrial variants

- **MitorSaw** -- mitochondrial variant caller and haplotype statistics.<br>
  [Documentation](https://github.com/PacificBiosciences/mitorsaw)<br>
  Citation: [PacificBiosciences. MitorSaw. GitHub repository.](https://github.com/PacificBiosciences/mitorsaw)

## Paralog genes

- **Paraphase** -- paralog-aware variant calling for genes with high sequence similarity (e.g., _SMN1_/_SMN2_, _PMS2_/_PMS2CL_).<br>
  [Documentation](https://github.com/PacificBiosciences/paraphase)<br>
  Citation: [Chen, X., et al. (2025). Genome-wide profiling of highly similar paralogous genes using HiFi sequencing. Nature Communications.](https://doi.org/10.1038/s41467-025-57505-2)<br>
  For the _SMN1_/_SMN2_ application specifically: [Chen, X., et al. (2023). Comprehensive SMN1 and SMN2 profiling for spinal muscular atrophy analysis using long-read PacBio HiFi sequencing. The American Journal of Human Genetics.](https://doi.org/10.1016/j.ajhg.2023.01.001)
- **Kivvi** -- targeted genotyping for the KIV2 (LPA) and D4Z4 (FSHD) repeat regions.<br>
  [Documentation](https://github.com/PacificBiosciences/kivvi)<br>
  Citation: [Chen, X., et al. (2026). Resolution of the D4Z4 repeat responsible for facioscapulohumeral muscular dystrophy with HiFi sequencing. bioRxiv, 2026.04.10.717730.](https://doi.org/10.64898/2026.04.10.717730) (pre-print, covers the D4Z4/FSHD side of Kivvi; no dedicated publication for the KIV2/LPA side)

## Pharmacogenomics

- **StarPhase** -- HLA typing and pharmacogenomic (PGx) diplotype calling.<br>
  [Documentation](https://github.com/PacificBiosciences/pb-StarPhase/blob/main/docs/user_guide.md)<br>
  Citation: [Holt, J. M., et al. (2024). StarPhase: comprehensive phase-aware pharmacogenomic diplotyper for long-read sequencing data. bioRxiv, 2024.12.10.627527.](https://doi.org/10.1101/2024.12.10.627527) (pre-print; no peer-reviewed publication yet)

## QC and statistics

- **pbjam** -- read-length, quality, and alignment-QC statistics and plotting (`bam_stats`).<br>
  [Documentation](https://github.com/PacificBiosciences/pbjam)<br>
  Citation: [PacificBiosciences. pbjam. GitHub repository.](https://github.com/PacificBiosciences/pbjam)

## Reference data

- **Reference data bundle** -- packaged reference FASTA, TRGT/Sawfish/MethBat region files, and per-reference defaults, distributed as [reference data containers](reference_data_containers.md) per `ref_name`.<br>
  [Documentation](reference_data_containers.md)<br>
  Citation: [Lake, J. A., Wenger, A. M., & Rowell, W. J. (2026). Reference data bundle for PacificBiosciences/HiFi-human-WGS-WDL (Version 4.0.0) [Data set]. Zenodo.](https://doi.org/10.5281/zenodo.21517827)

## Utility tools and libraries

These ship in the `pb_wdl_base` container (see [Tool versions and Containers](tools_containers.md)) and back many tasks across the pipeline (VCF/BAM manipulation, region arithmetic, and the Python plotting/QC scripts) rather than being scientific-analysis tools in their own right.

- **HTSlib** -- core C library for SAM/BAM/CRAM/VCF/BCF I/O, underlying samtools, bcftools, and pysam.<br>
  [Documentation](https://github.com/samtools/htslib)<br>
  Citation: [Bonfield, J. K., et al. (2021). HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience, 10(2), giab007.](https://doi.org/10.1093/gigascience/giab007)
- **samtools** / **bcftools** -- BAM/CRAM and VCF/BCF manipulation, respectively.<br>
  [samtools docs](http://www.htslib.org/doc/samtools.html) / [bcftools docs](http://www.htslib.org/doc/bcftools.html)<br>
  Citation: [Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008.](https://doi.org/10.1093/gigascience/giab008)
- **bedtools** -- genomic interval arithmetic (intersect, merge, etc.).<br>
  [Documentation](https://bedtools.readthedocs.io/)<br>
  Citation: [Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.](https://doi.org/10.1093/bioinformatics/btq033)
- **NumPy** -- array/numerical computing, used by several QC/plotting scripts.<br>
  [Documentation](https://numpy.org/doc/)<br>
  Citation: [Harris, C. R., et al. (2020). Array programming with NumPy. Nature, 585, 357-362.](https://doi.org/10.1038/s41586-020-2649-2)
- **pandas** -- tabular data handling in QC/reporting scripts.<br>
  [Documentation](https://pandas.pydata.org/docs/)<br>
  Citation: [McKinney, W. (2010). Data structures for statistical computing in Python. Proceedings of the 9th Python in Science Conference, 56-61.](https://doi.org/10.25080/Majora-92bf1922-00a)
- **Matplotlib** -- plotting backend for generated QC figures.<br>
  [Documentation](https://matplotlib.org/stable/)<br>
  Citation: [Hunter, J. D. (2007). Matplotlib: a 2D graphics environment. Computing in Science & Engineering, 9(3), 90-95.](https://doi.org/10.1109/MCSE.2007.55)
- **seaborn** -- statistical plotting on top of Matplotlib.<br>
  [Documentation](https://seaborn.pydata.org/)<br>
  Citation: [Waskom, M. L. (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021.](https://doi.org/10.21105/joss.03021)
- **pysam** -- Python bindings for HTSlib, used by several scripts to read BAM/VCF files.<br>
  [Documentation](https://pysam.readthedocs.io/)<br>
  Citation: [pysam-developers. pysam. GitHub repository.](https://github.com/pysam-developers/pysam) (pysam has no dedicated publication; its own FAQ recommends citing HTSlib and SAMtools/BCFtools above instead)
