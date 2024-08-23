eff=0.95
layer=6
noise=1
sbatch --output output/FWUV_eff${eff}_layer${layer}_noise${noise}.out FWUV.sh $eff $layer $noise

#sbatch --output output/Sig_FWUV_eff${eff}_layer${layer}_noise${noise}.out FWUV.sh $eff $layer $noise
