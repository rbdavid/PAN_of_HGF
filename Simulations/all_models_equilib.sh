#!/bin/bash
#BSUB -P BIP198
#BSUB -W 2:00
#BSUB -nnodes 1
#BSUB -alloc_flags gpudefault
#BSUB -J openmm
#BSUB -o hgf_wt.%J.out
#BSUB -e hgf_wt.%J.err

module load cuda/11.0.3 gcc/10.2.0

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

mkdir -p $MODEL_1_DIR/nvt_equilib
mkdir -p $MODEL_2_DIR/nvt_equilib
mkdir -p $MODEL_3_DIR/nvt_equilib
mkdir -p $MODEL_4_DIR/nvt_equilib
mkdir -p $MODEL_5_DIR/nvt_equilib

date
echo "RUNNING NVT EQUILIB"
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_nvt_simulations.py $SUBMIT_DIR/model_1/hgf_36_125_wt_81b40_relaxed_model_1_sim_ready.prmtop $SUBMIT_DIR/model_1/model_1_min.2/model_1.min.2.rst7 $MODEL_1_DIR/nvt_equilib/ > $SUBMIT_DIR/model_1/nvt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_nvt_simulations.py $SUBMIT_DIR/model_2/hgf_36_125_wt_81b40_relaxed_model_2_sim_ready.prmtop $SUBMIT_DIR/model_2/model_2_min.2/model_2.min.2.rst7 $MODEL_2_DIR/nvt_equilib/ > $SUBMIT_DIR/model_2/nvt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_nvt_simulations.py $SUBMIT_DIR/model_3/hgf_36_125_wt_81b40_relaxed_model_3_sim_ready.prmtop $SUBMIT_DIR/model_3/model_3_min.2/model_3.min.2.rst7 $MODEL_3_DIR/nvt_equilib/ > $SUBMIT_DIR/model_3/nvt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_nvt_simulations.py $SUBMIT_DIR/model_4/hgf_36_125_wt_81b40_relaxed_model_4_sim_ready.prmtop $SUBMIT_DIR/model_4/model_4_min.2/model_4.min.2.rst7 $MODEL_4_DIR/nvt_equilib/ > $SUBMIT_DIR/model_4/nvt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_nvt_simulations.py $SUBMIT_DIR/model_5/hgf_36_125_wt_81b40_relaxed_model_5_sim_ready.prmtop $SUBMIT_DIR/model_5/model_5_min.2/model_5.min.2.rst7 $MODEL_5_DIR/nvt_equilib/ > $SUBMIT_DIR/model_5/nvt_equilib.out &

wait

mkdir -p $MODEL_1_DIR/npt_equilib
mkdir -p $MODEL_2_DIR/npt_equilib
mkdir -p $MODEL_3_DIR/npt_equilib
mkdir -p $MODEL_4_DIR/npt_equilib
mkdir -p $MODEL_5_DIR/npt_equilib

echo "RUNNING NPT EQUILIB"
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_npt_simulations.py $SUBMIT_DIR/model_1/hgf_36_125_wt_81b40_relaxed_model_1_sim_ready.prmtop $MODEL_1_DIR/nvt_equilib/hgf_36_125_wt_81b40_relaxed_model_1_sim_ready.nvt_equilib.chkpt $MODEL_1_DIR/npt_equilib/ > $SUBMIT_DIR/model_1/npt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_npt_simulations.py $SUBMIT_DIR/model_2/hgf_36_125_wt_81b40_relaxed_model_2_sim_ready.prmtop $MODEL_2_DIR/nvt_equilib/hgf_36_125_wt_81b40_relaxed_model_2_sim_ready.nvt_equilib.chkpt $MODEL_2_DIR/npt_equilib/ > $SUBMIT_DIR/model_2/npt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_npt_simulations.py $SUBMIT_DIR/model_3/hgf_36_125_wt_81b40_relaxed_model_3_sim_ready.prmtop $MODEL_3_DIR/nvt_equilib/hgf_36_125_wt_81b40_relaxed_model_3_sim_ready.nvt_equilib.chkpt $MODEL_3_DIR/npt_equilib/ > $SUBMIT_DIR/model_3/npt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_npt_simulations.py $SUBMIT_DIR/model_4/hgf_36_125_wt_81b40_relaxed_model_4_sim_ready.prmtop $MODEL_4_DIR/nvt_equilib/hgf_36_125_wt_81b40_relaxed_model_4_sim_ready.nvt_equilib.chkpt $MODEL_4_DIR/npt_equilib/ > $SUBMIT_DIR/model_4/npt_equilib.out &
jsrun --smpiargs="none" -n 1 -a 1 -c 1 -g 1 python3 $SUBMIT_DIR/run_npt_simulations.py $SUBMIT_DIR/model_5/hgf_36_125_wt_81b40_relaxed_model_5_sim_ready.prmtop $MODEL_5_DIR/nvt_equilib/hgf_36_125_wt_81b40_relaxed_model_5_sim_ready.nvt_equilib.chkpt $MODEL_5_DIR/npt_equilib/ > $SUBMIT_DIR/model_5/npt_equilib.out &

wait

mv $MODEL_1_DIR/* $SUBMIT_DIR/model_1/
mv $MODEL_2_DIR/* $SUBMIT_DIR/model_2/
mv $MODEL_3_DIR/* $SUBMIT_DIR/model_3/
mv $MODEL_4_DIR/* $SUBMIT_DIR/model_4/
mv $MODEL_5_DIR/* $SUBMIT_DIR/model_5/

