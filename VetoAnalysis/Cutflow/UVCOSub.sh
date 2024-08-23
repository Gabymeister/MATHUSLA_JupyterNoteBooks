eff=0.95
layer=6
noise=1
sbatch --output output/UVCO_eff${eff}_layer${layer}_noise${noise}.out UVCO.sh $eff $layer $noise

#sbatch --output output/Sig_UVCO_eff${eff}_layer${layer}_noise${noise}.out UVCO.sh $eff $layer $noise
