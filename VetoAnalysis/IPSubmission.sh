effs=(0.90 0.95)
layers=(4 6)
noises=(1 10)
distances=(100 200 300 500 700 1000)
for eff in ${effs[@]}; do
	for layer in ${layers[@]}; do
		for noise in ${noises[@]}; do
			for dist in ${distances[@]}; do
				sbatch --output output/IP_eff${eff}_layer${layer}_noise${noise}_dist${dist}.out IPVeto.sh $eff $layer $noise ${dist}
			done
		done
	done
done

#for eff in ${effs[@]}; do
#	for layer in ${layers[@]}; do
#		for noise in ${noises[@]}; do
#			for dist in ${distances[@]}; do
#				sbatch --output output/Sig_IP_eff${eff}_layer${layer}_noise${noise}_dist${dist}.out IPVeto.sh $eff $layer $noise ${dist}
#			done
#		done
#	done
#done
