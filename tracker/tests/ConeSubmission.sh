
#effs=(0.90 0.95 1.00)
#layers=(4 5 6)
#noises=(0.1 1 10)

#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			sbatch --output output/Cone_eff${eff}_layer${layer}_noise${noise}.out ConeVeto.sh $eff $layer $noise
#		done
#	done
#done

effs=(0.90 0.95 1.00)
layers=(4)
noises=(0.1 1 10)

for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			sbatch --output output/Sig_Cone_eff${eff}_layer${layer}_noise${noise}.out ConeVeto.sh $eff $layer $noise
		done
	done
done
