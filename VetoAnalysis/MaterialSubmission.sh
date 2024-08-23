effs=(0.90 0.95)
layers=(4 6)
noises=(1 10)
distances=(1 5 10 25 50 100)
stricts=("False" "True")
for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			for dist in ${distances[@]}; do
				for strict in ${stricts}; do
					sbatch --output output/MAT_eff${eff}_layer${layer}_noise${noise}_dist${dist}_strict${strict}.out MaterialVeto.sh $eff $layer $noise $dist $strict
				done
			done
		done
	done
done

#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			for dist in ${distances[@]}; do
#				for strict in ${stricts}; do
#					sbatch --output output/Sig_MAT_eff${eff}_layer${layer}_noise${noise}_dist${dist}_strict${strict}.out MaterialVeto.sh $eff $layer $noise $dist $strict
#				done
#			done
#		done
#	done
#done
