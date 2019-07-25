#!/bin/bash
#SBATCH --job-name=nf-basicstats
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=control,normal
#################

nextflow run main.nf -resume -c nextflow.config --samples samples.tsv -with-timeline -with-trace -with-dag pipeline.png
