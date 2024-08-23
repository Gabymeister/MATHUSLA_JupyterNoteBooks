eff=0.95
layer=6
noise=1
sbatch --output output/TNCO_eff${eff}_layer${layer}_noise${noise}.out TNCO.sh $eff $layer $noise

#sbatch --output output/Sig_TNCO_eff${eff}_layer${layer}_noise${noise}.out TNCO.sh $eff $layer $noise
