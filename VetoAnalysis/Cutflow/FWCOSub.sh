eff=0.95
layer=6
noise=1
sbatch --output output/FWCO_eff${eff}_layer${layer}_noise${noise}.out FWCO.sh $eff $layer $noise

#sbatch --output output/Sig_FWCO_eff${eff}_layer${layer}_noise${noise}.out FWCO.sh $eff $layer $noise
