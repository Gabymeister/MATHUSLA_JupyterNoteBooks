eff=0.95
layer=6
noise=1
sbatch --output output/BestCase_FWIPUVMTCO_eff${eff}_layer${layer}_noise${noise}.out BestCase.sh $eff $layer $noise

#sbatch --output output/Sig_BestCase_FWIPUVMTCO_eff${eff}_layer${layer}_noise${noise}.out BestCase.sh $eff $layer $noise
