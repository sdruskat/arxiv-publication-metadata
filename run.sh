#!/bin/bash
echo "Running Snakemake workflow with 12 cores, keeping local storage."
snakemake --keep-storage-local-copies --software-deployment-method conda -c 12
echo "Done."

