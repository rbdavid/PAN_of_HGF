#!/bin/bash 

NMODELS=5
NPRODS=19

for ((model=1;model<=$NMODELS;model+=1))
do
	mkdir model_$model
	echo $model
	for ((prod=0;prod<=$NPRODS;prod+=5))
	do
		((a=$prod+4))
		echo $prod $a
		time python3 calc_avg_structure.py ../../truncated.pdb ../../model_$model/trajectories/ $model $prod $a average_list.00.$NPRODS.output
	done
	time python3 weighted_average_models.py average_list.00.$NPRODS.output ../../truncated.pdb
	mv *output model_$model/
	mv *avg_structure.pdb model_$model/
	mv *average_structure.pdb model_$model/
done

