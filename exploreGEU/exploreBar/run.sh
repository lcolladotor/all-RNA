#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=20G,h_fsize=4G
#$ -N exploreBar

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p /dcs01/ajaffe/Brain/derRuns/all-RNA/exploreGEU/exploreBar/logs

# Make plots
cd /dcs01/ajaffe/Brain/derRuns/all-RNA/exploreGEU/exploreBar
module load R/3.1.x
Rscript exploreTop_barcode.R

## Move log files into the logs directory
mv /dcs01/ajaffe/Brain/derRuns/all-RNA/exploreGEU/exploreBar/exploreBar.* /dcs01/ajaffe/Brain/derRuns/all-RNA/exploreGEU/exploreBar/logs/

echo "**** Job ends ****"
date