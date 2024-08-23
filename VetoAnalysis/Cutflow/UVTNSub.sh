eff=0.95
layer=6
noise=1
sbatch --output output/UVTN_eff${eff}_layer${layer}_noise${noise}.out UVTN.sh $eff $layer $noise

#sbatch --output output/Sig_UVTN_eff${eff}_layer${layer}_noise${noise}.out UVTN.sh $eff $layer $noise
