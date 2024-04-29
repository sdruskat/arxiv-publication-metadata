<!--
SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>

SPDX-License-Identifier: CC0-1.0
-->

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)  
[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)

# Snakemake workflow: Extract LUT from ArXiv OAI-PMH XML

Snakemake workflow to extract metadata from ArXiv OAI-PMH XML 
harvested with [`metha`](https://github.com/miku/metha),
and write it to JSON lookup tables for better accessibility.

## Documentation

The technical documentation and description of outputs is in [workflow/documentation.md](workflow/documentation.md).

## Running the workflow

You need to have `conda` installed to create and activate a new environment.

```bash
conda env create -n arxiv-metadata --file conda-environment.yaml
conda activate arxiv-metadata
```

You also need to get an access token from [Zenodo](https://zenodo.org) and set it to the following
two environment variables: 

```shell
export SNAKEMAKE_STORAGE_ZENODO_ACCESS_TOKEN=<Your Zenodo access token>
export SNAKEMAKE_STORAGE_ZENODO_RESTRICTED_ACCESS_TOKEN=<Your Zenodo access token>
```

Run with `-–keep-storage-local-copies` to avoid downloading resources over and over again.
Also run with `--software-deployment-method conda` to use global conda packages.

```shell
snakemake --keep-storage-local-copies --software-deployment-method conda -c <number-of-cores-to-use>
```

You can use [`run.sh`](run.sh) to run the workflow this way, and with 12 cores.
