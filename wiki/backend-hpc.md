# Installing and configuring for HPC backends

Either `miniwdl` or `Cromwell` can be used to run workflows on the HPC.

## Installing and configuring `miniwdl`

### Requirements

- [`miniwdl`](https://github.com/chanzuckerberg/miniwdl) >= 1.9.0
- [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm)

### Configuration

An [example miniwdl.cfg file](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/hpc/miniwdl.cfg) is provided here. This should be placed at `~/.config/miniwdl.cfg` and edited to match your slurm configuration. This allows running workflows using a basic SLURM setup.

## Installing and configuring `Cromwell`

Cromwell supports a number of different HPC backends; see [Cromwell's documentation](https://cromwell.readthedocs.io/en/stable/backends/HPC/) for more information on configuring each of the backends.  Cromwell can be used in a standalone "run" mode, or in "server" mode to allow for multiple users to submit workflows.  In the example below, we provide example commands for running Cromwell in "run" mode.

## Running the workflow

### Filling out workflow inputs

Fill out any information missing in [the inputs file](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/backends/hpc/singleton.hpc.inputs.json). Once you have downloaded the reference data bundle, ensure that you have replaced the `<local_path_prefix>` in the input template file with the local path to the reference datasets on your HPC.

See [the inputs section of the singleton README](./singleton#inputs) for more information on the structure of the inputs.json file.

#### Running via miniwdl

```bash
miniwdl run workflows/singleton.wdl --input <inputs_json_file>
```

#### Running via Cromwell

```bash
cromwell run workflows/singleton.wdl --input <inputs_json_file>
```

## Reference data bundle

[<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14027047.svg" alt="10.5281/zenodo.14027047">](https://zenodo.org/records/14027047)

Reference data is hosted on Zenodo at [10.5281/zenodo.14027047](https://zenodo.org/record/14027047). Download the reference data bundle and extract it to a location on your HPC, then update the input template file with the path to the reference data.

```bash
## download the reference data bundle
wget https://zenodo.org/record/14027047/files/hifi-wdl-resources-v2.0.0.tar

## extract the reference data bundle and rename as dataset
tar -xvf hifi-wdl-resources-v2.0.0.tar
```
