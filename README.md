[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)  
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

# Snakemake workflow: Extract LUT from ArXiv OAI-PMH XML

Snakemake workflow to extract metadata from ArXiv OAI-PMH XML 
harvested with [`metha`](https://github.com/miku/metha),
and write it to a JSON lookup table for better accessibility.

## Running the workflow

Run with `-–keep-storage-local-copies` to avoid downloading resources over and over again.

You can use [`run.sh`](run.sh) to run the workflow this way, and with 12 cores.
