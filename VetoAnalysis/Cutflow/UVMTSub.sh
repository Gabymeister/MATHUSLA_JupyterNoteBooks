eff=0.95
layer=6
noise=1
sbatch --output output/UVMT_eff${eff}_layer${layer}_noise${noise}.out UVMT.sh $eff $layer $noise

#sbatch --output output/Sig_UVMT_eff${eff}_layer${layer}_noise${noise}.out UVMT.sh $eff $layer $noise
