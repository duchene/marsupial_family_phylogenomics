#!/bin/bash
#PBS -P RDS-FSC-Phylogenomics-RW
#PBS -l select=1:ncpus=1:mem=16GB
#PBS -l walltime=24:00:00
#PBS -q defaultQ
cd $PBS_O_WORKDIR
module load R
module load phyml
Rscript run.mladeq.Rscript
