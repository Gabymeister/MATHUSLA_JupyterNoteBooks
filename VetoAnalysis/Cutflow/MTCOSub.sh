eff=0.95
layer=6
noise=1
sbatch --output output/MTCO_eff${eff}_layer${layer}_noise${noise}.out MTCO.sh $eff $layer $noise

#sbatch --output output/Sig_MTCO_eff${eff}_layer${layer}_noise${noise}.out MTCO.sh $eff $layer $noise
