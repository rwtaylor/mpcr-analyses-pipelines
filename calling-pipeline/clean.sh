#!/bin/bash
# Cleans up pipeline for a re-run.
# Removes all work in progrss - but does not delete final outputs.
rm -r work
rm -r .nextflow
rm .nextflow*
rm dag*
rm trace*
rm timeline*
rm slurm*
rm pipeline.png*
