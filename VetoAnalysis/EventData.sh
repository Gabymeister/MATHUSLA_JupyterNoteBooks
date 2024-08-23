#!/bin/bash
#SBATCH --time=5:30:00
#SBATCH --account=rrg-mdiamond
#SBATCH --mem=2G

module use --append /project/def-mdiamond/soft/modules

module load scdms/V05-00

echo running Event Data cutflow
singularity --silent exec --cleanenv -B /project/rrg-mdiamond/owhgabri/ -B /project/rrg-mdiamond/data/MATHUSLA/simulation/ /project/def-mdiamond/soft/releases/cdmsfull_V05-00.sif python3 EventData.py ${1} ${2} ${3}
