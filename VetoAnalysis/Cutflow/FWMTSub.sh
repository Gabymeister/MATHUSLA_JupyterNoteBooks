eff=0.95
layer=6
noise=1
sbatch --output output/FWMT_eff${eff}_layer${layer}_noise${noise}.out FWMT.sh $eff $layer $noise

#sbatch --output output/Sig_FWMT_eff${eff}_layer${layer}_noise${noise}.out FWMT.sh $eff $layer $noise
