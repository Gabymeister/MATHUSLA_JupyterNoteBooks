
effs=(0.90 0.95)
layers=(4 6)
noises=(1 10)
chi2s=(1 5 7 10 15 20 50 100)
distances=(100 200 300 500 700 1000)
for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			for chi2 in ${chi2s[@]}; do
				sbatch --output output/FW_eff${eff}_layer${layer}_noise${noise}_chi2${chi2}.out FW_Veto.sh $eff $layer $noise $chi2 0
			done
			for dist in ${distances[@]}; do
				sbatch --output output/FW_eff${eff}_layer${layer}_noise${noise}_dist${dist}.out FW_Veto.sh $eff $layer $noise 0 $dist
			done
		done
	done
done


#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			for chi2 in ${chi2s[@]}; do
#				sbatch --output output/Sig_FW_eff${eff}_layer${layer}_noise${noise}_chi2${chi2}.out FW_Veto.sh $eff $layer $noise $chi2 0
#			done
#			for dist in ${distances[@]}; do
#				sbatch --output output/Sig_FW_eff${eff}_layer${layer}_noise${noise}_dist${dist}.out FW_Veto.sh $eff $layer $noise 0 $dist
#			done
#		done
#	done
#done
