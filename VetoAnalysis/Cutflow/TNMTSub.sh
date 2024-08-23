eff=0.95
layer=6
noise=1
sbatch --output output/TNMT_eff${eff}_layer${layer}_noise${noise}.out TNMT.sh $eff $layer $noise

#sbatch --output output/Sig_TNMT_eff${eff}_layer${layer}_noise${noise}.out TNMT.sh $eff $layer $noise
