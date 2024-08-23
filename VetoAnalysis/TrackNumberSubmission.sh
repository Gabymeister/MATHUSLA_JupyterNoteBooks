effs=(0.90 0.95)
layers=(4 6)
noises=(1 10)
nums=(3 4)
stricts=("False" "True")
for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			for num in ${nums[@]}; do
				for strict in ${stricts[@]}; do
					sbatch --output output/TN_eff${eff}_layer${layer}_noise${noise}_TN${num}_strict${strict}.out TrackNumberVeto.sh $eff $layer $noise $num $strict
				done
			done
		done
	done
done

#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			for num in ${nums[@]}; do
#				for strict in ${stricts[@]}; do
#					sbatch --output output/Sig_TN_eff${eff}_layer${layer}_noise${noise}_TN${num}_strict${strict}.out TrackNumberVeto.sh $eff $layer $noise $num $strict
#				done
#			done
#		done
#	done
#done

