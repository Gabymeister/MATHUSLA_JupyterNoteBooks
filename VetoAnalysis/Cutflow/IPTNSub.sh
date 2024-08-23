eff=0.95
layer=6
noise=1
sbatch --output output/IPTN_eff${eff}_layer${layer}_noise${noise}.out IPTN.sh $eff $layer $noise

#sbatch --output output/Sig_IPTN_eff${eff}_layer${layer}_noise${noise}.out IPTN.sh $eff $layer $noise
