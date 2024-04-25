#!/bin/bash
echo "Running Snakemake workflow with 12 cores, keeping local storage."
snakemake --keep-storage-local-copies -c 12
echo "Done."

