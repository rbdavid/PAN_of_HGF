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
		time python3 calc_avg_structure.py /media/rbdavid/edi/Projects/ornl_HGF/hgf_36_125_wt/truncated.pdb /media/rbdavid/edi/Projects/ornl_HGF/hgf_36_125_wt/model_$model/trajectories/ $model $prod $a average_list.00.$NPRODS.output
	done
	time python3 weighted_average.py average_list.00.$NPRODS.output /media/rbdavid/edi/Projects/ornl_HGF/hgf_36_125_wt/truncated.pdb
	mv *output model_$model/
	mv *avg_structure.pdb model_$model/
	mv *average_structure.pdb model_$model/
done	
