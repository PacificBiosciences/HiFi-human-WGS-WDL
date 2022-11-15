# wdl-humanwgs

PacBio's human WGS workflow written in [Workflow Description Language (WDL)](https://openwdl.org/).

For the snakemake version of these workflows, see [here](https://github.com/PacificBiosciences/pb-human-wgs-workflow-snakemake).

Docker images used by these workflows are defined [here](https://github.com/PacificBiosciences/wdl-dockerfiles).


## Workflows

### Main

Calls all subworkflows for a full analysis.

**Workflow**: [workflows/main.wdl](workflows/main.wdl)
**Inputs**: [workflows/inputs.json](workflows/inputs.json)


### SMRT Cell analysis

Aligns reads to a reference genome and generates statistics on alignment depth, read length, and alignment quality.

**Workflow**: [workflows/smrtcell_analysis/smrtcell_analysis.wdl](workflows/smrtcell_analysis/smrtcell_analysis.wdl)
**Inputs**: [workflows/smrtcell_analysis/inputs.json](workflows/smrtcell_analysis/inputs.json)


### Sample analysis

Calls and phases small and large variants and structural variants.

**Workflow**: [workflows/sample_analysis/sample_analysis.wdl](workflows/sample_analysis/sample_analysis.wdl)
**Inputs**: [workflows/sample_analysis.inputs.json](workflows/sample_analysis/inputs.json)

### [Common](workflows/common)

These are resources that are used across workflows, including structs and common tasks.
