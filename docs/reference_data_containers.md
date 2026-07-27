# Reference data containers

Static workflow assets that used to be passed in as individual input files — the reference FASTA and index, TRGT/Sawfish/MethBat region files, and per-reference defaults like `paraphase_genome_build` — are now packaged together into a single container image per reference genome. The raw files backing these containers are archived at [Zenodo](https://zenodo.org/records/21517827).

## How it works

`singleton.wdl` and `family.wdl` each define a `reference_container` map keyed by `ref_name`:

| `ref_name` | Container |
| ---------- | --------- |
| `GRCh38` | [workflow-data-container-hifi-human-wgs-wdl-grch38@sha256:5e3f23e44c1c09762838e81677f7d136367be2a4063efc608129841cef745b42](https://quay.io/repository/pacbio/workflow-data-container-hifi-human-wgs-wdl-grch38/manifest/sha256:5e3f23e44c1c09762838e81677f7d136367be2a4063efc608129841cef745b42) (v4.0.0) |
| `GRCh38_GIABv3` | [workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3@sha256:f3799c1eff1b816a95ea06ce567d3324d06e945428567077d96eca07c76ad0aa](https://quay.io/repository/pacbio/workflow-data-container-hifi-human-wgs-wdl-grch38_giabv3/manifest/sha256:f3799c1eff1b816a95ea06ce567d3324d06e945428567077d96eca07c76ad0aa) (v4.0.0) |

Like our tool containers (see [Tool versions and Containers](tools_containers.md)), these are hosted on [quay.io/pacbio](https://quay.io/repository/pacbio) and referenced by their sha256 digest for reproducibility and call-caching. Adding support for a new reference means building and pushing a new container, then adding a `ref_name` entry to this map — no other WDL changes are needed.

The `unpack_container_manifest` task pulls the selected image and runs `/opt/scripts/unpack_container.py` against a `manifest.json` baked into the image at `/opt/manifests/manifest.json`, which describes and locates all of the container's static workflow inputs. This task is marked `cacheable: true`, an engine-specific hint (currently honored only by `sprocket`) so the unpack only has to happen once per reference rather than once per run; other engines will simply re-run it each time.

## What's inside

Unpacking a reference container produces:

- `ref_fasta`, `ref_index` — the reference genome FASTA and its index
- `max_norm_female_chrY_depth` — maximum expected normalized chrY depth for samples without chrY
- `trgt_tandem_repeat_bed` — TRGT tandem repeat catalog
- `sawfish_exclude_bed`(`_index`) — regions excluded from Sawfish CNV calling
- `sawfish_expected_bed_male`, `sawfish_expected_bed_female` — expected allosome copy number BEDs
- `methbat_region_tsv` — MethBat methylation profiling regions
- `paraphase_genome_build` — genome build parameter passed to Paraphase
- `run_starphase` — whether to run the StarPhase task for this reference
- `manifest_json` — the manifest itself, for provenance

Two of these can still be overridden per-run with a corresponding workflow input: `trgt_tandem_repeat_bed_override` and `methbat_region_tsv_override`. Everything else is fixed by the selected reference container.
