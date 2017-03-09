#!/bin/bash
#SBATCH --job-name=vcf-filter
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=hns,dpetrov
#################

nextflow run main.nf -resume -c nextflow.config -with-timeline -with-trace -with-dag dag.png
