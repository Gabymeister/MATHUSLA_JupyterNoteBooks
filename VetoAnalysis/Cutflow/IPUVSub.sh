eff=0.95
layer=6
noise=1
sbatch --output output/IPUV_eff${eff}_layer${layer}_noise${noise}.out IPUV.sh $eff $layer $noise

#sbatch --output output/Sig_IPUV_eff${eff}_layer${layer}_noise${noise}.out IPUV.sh $eff $layer $noise
