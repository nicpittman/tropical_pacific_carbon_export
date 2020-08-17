#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -l mem=8GB
#PBS -l ncpus=6
#PBS -N combine_moorings
#PBS -l storage=gdata/hh5+gdata/v45
#PBS -P v45
#PBS -q normal
#PBS -l wd 

module purge
module use /g/data3/hh5/public/modules/
module load conda/analysis3

python3 /g/data/v45/np1383/cb2/tropical_pacific_carbon_export/7b-gadi-Data_Combine_CO2_Physics_parallel.py 
