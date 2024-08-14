#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --mem=2G

module use --append /project/def-mdiamond/soft/modules

module load scdms/V05-00

echo running Cutoff veto 
singularity --silent exec --cleanenv -B /project/rrg-mdiamond/owhgabri/ /project/def-mdiamond/soft/releases/cdmsfull_V05-01.sif python3 CutoffVeto.py ${1} ${2} ${3} ${4}
