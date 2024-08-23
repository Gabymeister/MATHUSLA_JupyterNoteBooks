eff=0.95
layer=6
noise=1
sbatch --output output/FWTN_eff${eff}_layer${layer}_noise${noise}.out FWTN.sh $eff $layer $noise

#sbatch --output output/Sig_FWTN_eff${eff}_layer${layer}_noise${noise}.out FWTN.sh $eff $layer $noise
