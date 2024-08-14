
#effs=(0.90 0.95 1.00)
#layers=(4 5 6)
#noises=(0.1 1 10)
#stricts=("False" "True")
#
#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			for strict in ${stricts[@]}; do
#				sbatch --output output/Cutoff_eff${eff}_layer${layer}_noise${noise}_strict${strict}.out CutoffVeto.sh $eff $layer $noise $strict
#			done
#		done
#	done
#done

effs=(0.90 0.95 1.00)
layers=(4)
noises=(0.1 1 10)
stricts=("False" "True")

for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			for strict in ${stricts[@]}; do
				sbatch --output output/Sig_Cutoff_eff${eff}_layer${layer}_noise${noise}_strict${strict}.out CutoffVeto.sh $eff $layer $noise $strict
			done
		done
	done
done


