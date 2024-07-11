# Tool versions and Containers

Containers are used to package tools and their dependencies. This ensures that the tools are reproducible and can be run on any system that supports the container runtime.  Our containers are built using [Docker](https://www.docker.com/) and are compatible with any container runtime that supports the OCI Image Specification, like [Singularity](https://sylabs.io/singularity/) or [Podman](https://podman.io/).

Most of our containers are built on the `pb_wdl_base` container, which includes common bioinformatics tools and libraries.  We tag our containers with a version number and build count, but the containers are referenced within the WDL files by the sha256 sum tags for reproducibility and better compatibility with Cromwell and miniwdl call caching.

Our Dockerfiles can be inspected on GitHub, and the containers can be pulled from our [Quay.io organization](https://quay.io/pacbio).

We directly use `deepvariant`, `deepvariant-gpu`, `pharmcat`, and `glnexus` containers from their respective authors, although we have mirrored some for better compatibility with Cromwell call caching.

| Container | Major tool versions | Dockerfile | Container<br/>(links directly to `docker://` URI) |
| --------: | ------------------- | :---: | :---: |
| pb_wdl_base | <ul><li>htslib 1.20</li><li>bcftools 1.20</li><li>samtools 1.20</li><li>bedtools 2.31.0</li><li>python3.9</li><li>numpy 1.24.24</li><li>pandas 2.0.3</li><li>matplotlib 3.7.5</li><li>seaborn 0.13.2</li><li>pysam 0.22.1</li><li>vcfpy 0.13.8</li><li>biopython 1.83</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/pb_wdl_base) | [Container](docker://quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87) |
| pbmm2 | <ul><li>pbmm2 1.13.1</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/pbmm2) | [Container](docker://quay.io/pacbio/pbmm2@sha256:265eef770980d93b849d1ddb4a61ac449f15d96981054e91d29da89943084e0e) |
| mosdepth | <ul><li>mosdepth 0.3.8</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/mosdepth) | [Container](docker://quay.io/pacbio/mosdepth@sha256:f715c11100e9bb3562cce1c5e23a185cfcc92a6fec412b16c30c0250496cc0d1) |
| pbsv | <ul><li>pbsv 2.9.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/pbsv) | [Container](docker://quay.io/pacbio/pbsv@sha256:7626286e07dd185ca698efc80bd0d26cd3a139fe19781dfde5b6d07e895673cd) |
| trgt | <ul><li>trgt 1.0.0</li><li>`/opt/scripts/check_trgt_coverage.py` 0.1.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/1cbc3fb2d8ef1303b41a6d30cc2012d0de96fb10/docker/trgt) | [Container](docker://quay.io/pacbio/trgt@sha256:e626b0a102c11ad9d52bed5a5573052bc76560f3f02146c48babc4a76bc74f52) |
| hiphase | <ul><li>hiphase 1.4.2</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/baf058daa1f90435a6de9d35f6dff13169411870/docker/hiphase) | [Container](docker://quay.io/pacbio/hiphase@sha256:b366a0ddbe42b79e238941ce687b0e34aac8d8c1de9c2c67af18e53fcb6f0c69) |
| hificnv | <ul><li>hificnv 1.0.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/hificnv) | [Container](docker://quay.io/pacbio/hificnv@sha256:c9e2d07240299cfff655ae9a96eb604934879128bd7aed9e60af6619f6c36b9a) |
| paraphase | <ul><li>paraphase 3.1.1</li><li>minimap 2.28</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/paraphase) | [Container](docker://quay.io/pacbio/paraphase@sha256:a114ac5b9a682d7dc0fdf25c92cfb36f80c07ab4f1fb76b2e58092521b123a4d) |
| pbstarphase | <ul><li>pbstarphase 0.10.0</li><li>CPIC v0.9.0, 2024_04_04</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/pbstarphase) | [Container](docker://quay.io/pacbio/pbstarphase@sha256:578f99e86e977436420187bfa0ad820067e574cbfb00a5461c773c12e1ba29ad) |
| pb-cpg-tools | <ul><li>pb-cpg-tools 2.3.2</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/pb-cpg-tools) | [Container](docker://quay.io/pacbio/pb-cpg-tools@sha256:d6e63fe3f6855cfe60f573de1ca85fab27f4a68e24a7f5691a7a805a22af292d) |
| wgs_tertiary | <ul><li>`/opt/scripts/calculate_phrank.py` 2.0.0</li><li>`/opt/scripts/json2ped.py` 0.1.0</li></ul>Last built 2021-09-17:<ul><li>ensembl -> HGNC</li><li>ensembl -> HPO</li><li>HGNC -> inheritance</li><li>HPO DAG</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/wgs_tertiary) | [Container](docker://quay.io/pacbio/wgs_tertiary@sha256:8fc134fdf0665e14a67bf7a8b4b63f5ae891a370a1d50c9eec2059702440a3e2) |
| slivar | <ul><li>slivar 0.3.0</li><li>`/opt/scripts/add_comphet_phase.py` 0.1.0</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6b13cc246dd44e41903d17a660bb5432cdd18dbe/docker/slivar) | [Container](docker://quay.io/pacbio/slivar@sha256:35be557730d3ac9e883f1c2010fb24ac02631922f9b4948b0608d3e643a46e8b) |
| svpack | <ul><li>svpack 54b54db</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/6fc750b0c65b4a5c1eb65791eab9eed89864d858/docker/svpack) | [Container](docker://quay.io/pacbio/svpack@sha256:628e9851e425ed8044a907d33de04043d1ef02d4d2b2667cf2e9a389bb011eba) |
| deepvariant | <ul><li>DeepVariant 1.6.1</li></ul> |  | [Container](docker://google/deepvariant:1.6.1) |
| deepvariant-gpu | <ul><li>DeepVariant 1.6.1</li></ul> |  | [Container](docker://google/deepvariant:1.6.1-gpu) |
| pharmcat | <ul><li>PharmCat 2.12.0</li></ul> |  | [Container](docker://quay.io/pacbio/glnexus@sha256:ce6fecf59dddc6089a8100b31c29c1e6ed50a0cf123da9f2bc589ee4b0c69c8e) |
| glnexus | <ul><li>GLnexus 1.4.3</li></ul> |  | [Container](docker://pgkb/pharmcat:2.12.0) |