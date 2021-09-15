#!/bin/bash
#BSUB -P BIP198
#BSUB -W 2:00
#BSUB -nnodes 1
#BSUB -alloc_flags gpudefault
#BSUB -J openmm
#BSUB -o hgf_wt.%J.out
#BSUB -e hgf_wt.%J.err

module load cuda/11.0.3 gcc/10.2.0

prev_prod="$((aaa-1))"
prev_prod_zero=$(printf "%02d" $prev_prod)

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/ccs/home/davidsonrb/Apps/miniconda_SUMMIT/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/ccs/home/davidsonrb/Apps/miniconda_SUMMIT/etc/profile.d/conda.sh" ]; then
        . "/ccs/home/davidsonrb/Apps/miniconda_SUMMIT/etc/profile.d/conda.sh"
    else
        export PATH="/ccs/home/davidsonrb/Apps/miniconda_SUMMIT/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate openmm

export SUBMIT_DIR="/gpfs/alpine/bip198/proj-shared/HGF/WT"

echo "PREPPING MEMBERWORKS"

export MODEL_1_DIR="/gpfs/alpine/scratch/davidsonrb/bip198/HGF/WT/model_1"
export MODEL_2_DIR="/gpfs/alpine/scratch/davidsonrb/bip198/HGF/WT/model_2"
export MODEL_3_DIR="/gpfs/alpine/scratch/davidsonrb/bip198/HGF/WT/model_3"
export MODEL_4_DIR="/gpfs/alpine/scratch/davidsonrb/bip198/HGF/WT/model_4"
export MODEL_5_DIR="/gpfs/alpine/scratch/davidsonrb/bip198/HGF/WT/model_5"

date

mkdir -p $MODEL_1_DIR/prod_aaa
mkdir -p $MODEL_2_DIR/prod_aaa
mkdir -p $MODEL_3_DIR/prod_aaa
mkdir -p $MODEL_4_DIR/prod_aaa
mkdir -p $MODEL_5_DIR/prod_aaa

echo "RUNNING PRODUCTION RUN aaa"
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_production_sims.py $SUBMIT_DIR/model_1/hgf_36_125_wt_81b40_relaxed_model_1_sim_ready.prmtop $SUBMIT_DIR/model_1/prod_$prev_prod_zero/hgf_36_125_wt_81b40_relaxed_model_1_sim_ready.prod_$prev_prod_zero.chkpt $MODEL_1_DIR/prod_aaa/ aaa > $SUBMIT_DIR/model_1/prod_aaa.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_production_sims.py $SUBMIT_DIR/model_2/hgf_36_125_wt_81b40_relaxed_model_2_sim_ready.prmtop $SUBMIT_DIR/model_2/prod_$prev_prod_zero/hgf_36_125_wt_81b40_relaxed_model_2_sim_ready.prod_$prev_prod_zero.chkpt $MODEL_2_DIR/prod_aaa/ aaa > $SUBMIT_DIR/model_2/prod_aaa.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_production_sims.py $SUBMIT_DIR/model_3/hgf_36_125_wt_81b40_relaxed_model_3_sim_ready.prmtop $SUBMIT_DIR/model_3/prod_$prev_prod_zero/hgf_36_125_wt_81b40_relaxed_model_3_sim_ready.prod_$prev_prod_zero.chkpt $MODEL_3_DIR/prod_aaa/ aaa > $SUBMIT_DIR/model_3/prod_aaa.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_production_sims.py $SUBMIT_DIR/model_4/hgf_36_125_wt_81b40_relaxed_model_4_sim_ready.prmtop $SUBMIT_DIR/model_4/prod_$prev_prod_zero/hgf_36_125_wt_81b40_relaxed_model_4_sim_ready.prod_$prev_prod_zero.chkpt $MODEL_4_DIR/prod_aaa/ aaa > $SUBMIT_DIR/model_4/prod_aaa.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_production_sims.py $SUBMIT_DIR/model_5/hgf_36_125_wt_81b40_relaxed_model_5_sim_ready.prmtop $SUBMIT_DIR/model_5/prod_$prev_prod_zero/hgf_36_125_wt_81b40_relaxed_model_5_sim_ready.prod_$prev_prod_zero.chkpt $MODEL_5_DIR/prod_aaa/ aaa > $SUBMIT_DIR/model_5/prod_aaa.out &

wait

mv $MODEL_1_DIR/* $SUBMIT_DIR/model_1/
mv $MODEL_2_DIR/* $SUBMIT_DIR/model_2/
mv $MODEL_3_DIR/* $SUBMIT_DIR/model_3/
mv $MODEL_4_DIR/* $SUBMIT_DIR/model_4/
mv $MODEL_5_DIR/* $SUBMIT_DIR/model_5/

