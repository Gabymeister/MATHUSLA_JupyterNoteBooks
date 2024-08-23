effs=(0.90 0.95)
layers=(4 6)
noises=(1 10)
for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			sbatch --output output/Cone_eff${eff}_layer${layer}_noise${noise}.out ConeVeto.sh $eff $layer $noise
		done
	done
done

#effs=(0.90 0.95)
#layers=(4 6)
#noises=(1 10)
#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			sbatch --output output/Sig_Cone_eff${eff}_layer${layer}_noise${noise}.out ConeVeto.sh $eff $layer $noise
#		done
#	done
#done
