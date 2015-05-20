#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=5G,h_vmem=15G,h_fsize=10G
#$ -N testMean20
#$ -pe local 4
#$ -q shared.q@compute-061

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p logs

# Main script
R-bioc-devel -e "source('testMean.R')"

## Move log files into the logs directory
mv testMean20.* logs/

echo '**** Job ends ****'
date
