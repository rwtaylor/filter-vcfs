#!/bin/bash
#SBATCH --job-name=vcf-rsync
#SBATCH --time=2-0:0:0
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --partition=hns,dpetrov
#################

rsync -avh /scratch/PI/dpetrov/2016-tiger-wgs/vcfs/subsampled ryan@bio-dap15.stanford.edu:/zstor/2016-tiger-wgs/vcfs/
rsync -avh /scratch/PI/dpetrov/2016-tiger-wgs/vcfs/filtered ryan@bio-dap15.stanford.edu:/zstor/2016-tiger-wgs/vcfs/
rsync -avh /scratch/PI/dpetrov/2016-tiger-wgs/vcfs/raw ryan@bio-dap15.stanford.edu:/zstor/2016-tiger-wgs/vcfs/