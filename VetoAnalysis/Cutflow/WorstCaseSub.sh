eff=0.90
layer=4
noise=10
sbatch --output output/WorstCase_FWIPUVMTCO_eff${eff}_layer${layer}_noise${noise}.out WorstCase.sh $eff $layer $noise

#sbatch --output output/Sig_WorstCase_FWIPUVMTCO_eff${eff}_layer${layer}_noise${noise}.out WorstCase.sh $eff $layer $noise
