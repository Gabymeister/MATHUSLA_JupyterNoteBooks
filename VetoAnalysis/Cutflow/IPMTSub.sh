eff=0.95
layer=6
noise=1
sbatch --output output/IPMT_eff${eff}_layer${layer}_noise${noise}.out IPMT.sh $eff $layer $noise

#sbatch --output output/Sig_IPMT_eff${eff}_layer${layer}_noise${noise}.out IPMT.sh $eff $layer $noise
