eff=0.95
layer=6
noise=1
sbatch --output output/IPCO_eff${eff}_layer${layer}_noise${noise}.out IPCO.sh $eff $layer $noise

#sbatch --output output/Sig_IPCO_eff${eff}_layer${layer}_noise${noise}.out IPCO.sh $eff $layer $noise
