#!/bin/bash
#SBATCH -A bip198
#SBATCH -p batch
#SBATCH -t 6:00:00
#SBATCH -N 1
##SBATCH --mail-type=ALL
##SBATCH --mail-user=davidsonrb@ornl.gov
#SBATCH -J OFA
#SBATCH -o ./post_processing.out
#SBATCH -e ./post_processing.err

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ccs/home/davidsonrb/Apps/miniconda3_ANDES/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ccs/home/davidsonrb/Apps/miniconda3_ANDES/etc/profile.d/conda.sh" ]; then
        . "/ccs/home/davidsonrb/Apps/miniconda3_ANDES/etc/profile.d/conda.sh"
    else
        export PATH="/ccs/home/davidsonrb/Apps/miniconda3_ANDES/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate OpenFold-amber

DIRS=*/
#DIRS=model_1/
echo $DIRS
for dir in $DIRS
do
	echo $dir
	model_num=${dir:6:-1}

	cd $dir
	PRODS=prod_*/
	echo $PRODS
	
	mkdir truncated_trajs/
	
	sed -e s?AAA?model_"$model_num"_min.1/model_"$model_num".min.1.nc?g -e s?CCC?truncated_trajs/hgf_wt_model_"$model_num".min.1.trunc.dcd?g < ../cpptraj_truncate.inp > temp_truncate.inp
	cpptraj -p hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.prmtop -i temp_truncate.inp > truncated_trajs/cpptraj.min.1.log &
	sleep 5s
	
	sed -e s?AAA?model_"$model_num"_min.2/model_"$model_num".min.2.nc?g -e s?CCC?truncated_trajs/hgf_wt_model_"$model_num".min.2.trunc.dcd?g < ../cpptraj_truncate.inp > temp_truncate.inp
	cpptraj -p hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.prmtop -i temp_truncate.inp > truncated_trajs/cpptraj.min.2.log &
	sleep 5s
	
	sed -e s?AAA?nvt_equilib/hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.nvt_equilib.dcd?g -e s?CCC?truncated_trajs/hgf_wt_model_"$model_num".nvt_equilib.trunc.dcd?g < ../cpptraj_truncate.inp > temp_truncate.inp
	cpptraj -p hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.prmtop -i temp_truncate.inp > truncated_trajs/cpptraj.nvt_equilib.log &
	sleep 5s
	
	sed -e s?AAA?npt_equilib/hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.npt_equilib.dcd?g -e s?CCC?truncated_trajs/hgf_wt_model_"$model_num".npt_equilib.trunc.dcd?g < ../cpptraj_truncate.inp > temp_truncate.inp
	cpptraj -p hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.prmtop -i temp_truncate.inp > truncated_trajs/cpptraj.npt_equilib.log &
	sleep 5s
	
	for prod in $PRODS
	do
		echo $prod
		prod_num=${prod:5:-1}
		echo $prod_num
		sed -e s?AAA?prod_$prod_num/hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.prod_$prod_num.dcd?g -e s?CCC?truncated_trajs/hgf_wt_model_"$model_num".prod_$prod_num.trunc.dcd?g < ../cpptraj_truncate.inp > temp_truncate.inp
		cpptraj -p hgf_36_125_wt_81b40_relaxed_model_"$model_num"_sim_ready.prmtop -i temp_truncate.inp > truncated_trajs/cpptraj.prod_"$prod_num".log &
		sleep 5s
	done
	
	wait
	echo Done with truncating $dir files
	
	cd ../
done

