# wdl-humanwgs

PacBio's human WGS workflow written in [Workflow Description Language (WDL)](https://openwdl.org/).

For the snakemake version of these workflows, see [here](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake).

Docker images used by these workflows are defined [here](https://github.com/PacificBiosciences/wdl-dockerfiles).


## Workflows

### SMRT Cell analysis

Aligns reads to a reference genome and generates statistics on alignment depth, read length, and alignment quality.

**Workflow**: [workflows/smrtcell_analysis/smrtcell_analysis.wdl](workflows/smrtcell_analysis.wdl)
**Inputs**: [workflows/smrtcell_analysis/inputs.json]
