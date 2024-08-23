eff=0.95
layer=6
noise=1
sbatch --output output/FWIP_eff${eff}_layer${layer}_noise${noise}.out FWIP.sh $eff $layer $noise

#sbatch --output output/Sig_FWIP_eff${eff}_layer${layer}_noise${noise}.out FWIP.sh $eff $layer $noise
