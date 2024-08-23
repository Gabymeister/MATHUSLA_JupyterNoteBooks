#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --mem=2G

module use --append /project/def-mdiamond/soft/modules

module load scdms/V05-00

echo running MTCO
singularity --silent exec --cleanenv -B /project/rrg-mdiamond/owhgabri/ -B /project/rrg-mdiamond/data/MATHUSLA/simulation/ /project/def-mdiamond/soft/releases/cdmsfull_V05-00.sif python3 MTCO.py ${1} ${2} ${3}
