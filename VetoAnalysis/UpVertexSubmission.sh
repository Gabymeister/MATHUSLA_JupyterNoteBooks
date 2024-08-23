effs=(0.90 0.95)
layers=(4 6)
noises=(1 10)
for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			sbatch --output output/UpVertex_eff${eff}_layer${layer}_noise${noise}_chi2${chi2}.out UpVertexVeto.sh $eff $layer $noise
		done
	done
done

#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			sbatch --output output/Sig_UpVertex_eff${eff}_layer${layer}_noise${noise}_chi2${chi2}.out UpVertexVeto.sh $eff $layer $noise
#		done
#	done
#done

