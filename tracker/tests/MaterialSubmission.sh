effs=(0.90 0.95 1.00)
layers=(4 5 6)
noises=(0.1 1 10)
distances=(5 25 50)
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

#effs=(0.90 0.95 1.00)
#layers=(4)
#noises=(0.1 1 10)
#distances=(5 25 50)
#stricts=("False" "True")

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
